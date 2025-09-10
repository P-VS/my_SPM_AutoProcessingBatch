function params = my_spmbatch_ica1stlevel(sub,ses,run,task,datpath,params)

curdir = pwd;

matlabbatch = {};

%% Search for the data folders
ppparams = my_spmbatch_1stlevel_FindData(sub,ses,run,task,datpath,params);

if params.func.mruns && contains(params.func.use_runs,'separately')
    ppparams.resultmap = fullfile(ppparams.subpath,['ICA-' task '_' params.analysisname '_run-' num2str(run)]);
else
    ppparams.resultmap = fullfile(ppparams.subpath,['ICA-' task '_' params.analysisname]);
end

if params.do_ica && exist(ppparams.resultmap,'dir'); rmdir(ppparams.resultmap,'s'); end
if ~exist(ppparams.resultmap,'dir'), mkdir(ppparams.resultmap); end

Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(1).func(1).funcfile));
nvols = min([numel(Vfunc),50]);
fdata = spm_read_vols(Vfunc(1:nvols));
mask = my_spmbatch_mask(fdata);

Vmask = Vfunc(1);
rmfield(Vmask,'pinfo');
Vmask.fname = fullfile(ppparams.resultmap,['mask_' ppparams.frun(1).func(1).funcfile]);
Vmask.descrip = 'my_spmbatch - mask';
Vmask.dt = [spm_type('float32'),spm_platform('bigend')];
Vmask.n = [1 1];
Vmask = spm_write_vol(Vmask,mask);

ppparams.mask_file = Vmask.fname;

%% fMRI model specification
if params.tempsort_SPM_model
    tmp_resultmap = ppparams.resultmap;
    ppparams.resultmap = fullfile(tmp_resultmap,'SPMMAT');
    if exist(ppparams.resultmap,'dir'); rmdir(ppparams.resultmap,'s'); end
    mkdir(ppparams.resultmap)
    
    [matlabbatch,ppparams] = my_spmbatch_1stlevel_DefineModel(sub,ses,run,task,datpath,params,ppparams,matlabbatch);
    
    spm_jobman('run', matlabbatch)
    
    SPM_file = fullfile(ppparams.resultmap,'SPM.mat');
    
    %% Optimize GLM with TEDM
    if params.optimize_HRF
        SPM_file = my_spmmbatch_tedm(SPM_file,ppparams.resultmap,ppparams.mask_file);
        
        load(SPM_file)
    
        SPM.xX.X = SPM.TEDM.Param.Del;
    
        save(SPM_file,'SPM')
    end
    
    ppparams.resultmap = tmp_resultmap;
end

%% Adapt template parameter to specific subject data
% Get t_r
jsondat = jsondecode(fileread(ppparams.frun(1).funcjsonfile));
t_r = jsondat.("RepetitionTime");

for ir=1:numel(params.iruns)
    fsplit = split(ppparams.ppfmridat{ir}.sess{1}.func{1,1},'.nii');
    ica_source_file{ir} = [fsplit{1} '.nii'];
end

icatb_defaults;
icamat_file = fullfile(ppparams.resultmap,'ica__ica_parameter_info.mat');
compFiles = fullfile(ppparams.resultmap,'ica__sub01_component_ica_s1_.nii');

if params.do_ica || ~exist("icamat_file","file") || ~exist("compFiles","file")
    % Load template parameters .mat file
    icaparms_file = ica_1stlevel_make_gift_parameters(ica_source_file,t_r,numel(params.iruns),params,ppparams);
    
    %% Set up ICA
    icaparam_file = icatb_setup_analysis(icaparms_file);

    load(icaparam_file);

    %% Run Analysis (All steps)
    sesInfo = icatb_runAnalysis(sesInfo, [2:6]);
else
    load(icamat_file)
end

%%%%%%%%%% Get the required variables from sesInfo structure %%%%%%%%%%
% Number of subjects
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;

% Number of components
numComp = sesInfo.numComp;
dataType = sesInfo.dataType;
mask_ind = sesInfo.mask_ind;

% First scan
structFile = deblank(sesInfo.inputFiles(1).name(1, :));

%%%%%%%% End for getting the required vars from sesInfo %%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Get component data %%%%%%%%%%%%%%%%%

% Get the ICA Output files
icaOutputFiles = sesInfo.icaOutputFiles;

[subjectICAFiles, meanICAFiles, tmapICAFiles, meanALL_ICAFile] = icatb_parseOutputFiles('icaOutputFiles', icaOutputFiles, 'numOfSub', ...
        numOfSub, 'numOfSess', numOfSess, 'flagTimePoints', sesInfo.flagTimePoints);

% component files
if ~exist(meanALL_ICAFile.name)
    compFiles = subjectICAFiles(1).ses(1).name;
else
    compFiles = meanALL_ICAFile.name;
end

compFiles = icatb_fullFile('directory', ppparams.resultmap, 'files', compFiles);

utcompData = spm_read_vols(spm_vol(compFiles));
dim = size(utcompData);
HInfo.V = spm_vol(structFile);
HInfo.DIM = dim(1:3);
HInfo.VOX = double(HInfo.V(1).private.hdr.pixdim(2:4)); HInfo.VOX = abs(HInfo.VOX);

thcompData = permute(utcompData, [4 1 2 3]);
thcompData = reshape(thcompData, size(thcompData, 1), prod(dim(1:3)));
[thcompData] = icatb_applyDispParameters(thcompData, 1, 3, 1, dim(1:3), HInfo);
thcompData = reshape(thcompData, [size(thcompData,1), dim(1), dim(2), dim(3)]);
thcompData = permute(thcompData, [2 3 4 1]);

VC = spm_vol(compFiles);
for ic=1:numComp
    VC(ic).fname = fullfile(ppparams.resultmap, 'thres_compdata.nii'); 
    VC(ic).n = [ic 1];
    spm_write_vol(VC(ic), thcompData(:,:,:,ic));
end

compData = permute(utcompData, [4 1 2 3]);
compData = reshape(compData, size(compData, 1), prod(dim(1:3)));
[compData] = icatb_applyDispParameters(compData, 1, 1, 1, dim(1:3), HInfo);
compData = reshape(compData, [size(compData,1), dim(1), dim(2), dim(3)]);
compData = permute(compData, [2 3 4 1]);

% load time course
icaTimecourse = icatb_loadICATimeCourse(compFiles, 'real', [], [1:numComp]);

% Structural volume
dim = HInfo.DIM;
tdim = size(icaTimecourse(:,1));

%% Load model timecourses from design matrix
if params.tempsort_SPM_model
    load(SPM_file)
    
    modelTimecourse = SPM.xX.X;
    szX = size(modelTimecourse);
    
    useComp = zeros(szX(2),numComp);
    
    % Test correlation with model time course
    for ireg=1:szX(2)
        nsplit=split(SPM.xX.name{ireg},' ');
        if ~isempty(nsplit{2}) && ~contains(nsplit{2},'constant')
            regmap = zeros(HInfo.DIM(1:3));
            thregmap = zeros(HInfo.DIM(1:3));

            for ic=1:numComp
                p=polyfit(modelTimecourse(:,ireg),icaTimecourse(:,ic),1);
                f=polyval(p,modelTimecourse(:,ireg));
                SSE=sum((icaTimecourse(:,ic)-f).^2,'all');
                SST=sum((icaTimecourse(:,ic)-mean(icaTimecourse(:,ic),'all')).^2,'all');
                Rsqr(ireg,ic)=1-(SSE/SST);
    
                if Rsqr>0.4
                    useComp(ireg,ic) = 1; 
                    regmap = regmap+utcompData(:,:,:,ic); 
                    thregmap = thregmap+compData(:,:,:,ic);
                end
            end

            V = spm_vol(compFiles);
    
            Vout = V(1);
            Vout.fname = fullfile(ppparams.resultmap,['ica_map_' nsplit{2} '.nii']);
            Vout = spm_write_vol(Vout,regmap);
        
            Vout = V(1);
            Vout.fname = fullfile(ppparams.resultmap,['thres_ica_map_' nsplit{2} '.nii']);
            Vout = spm_write_vol(Vout,thregmap);
        
            clear Vout regmap thregmap
        end
    end

    save("components_per_SPMcondition.mat",'useComp','Rsqr')
end

%% Test overlap with atlas mask
if ~isempty(params.spatial_networks_masks)

    % Make brain and no-brain masks
    if ppparams.frun(1).func(1).funcfile(1)=='s', segfunc = ppparams.frun(1).func(1).funcfile(2:end); else segfunc = ppparams.frun(1).func(1).funcfile; end

    if ~isfile(fullfile(ppparams.preprocfmridir,['c1' segfunc])) || ~isfile(fullfile(ppparams.preprocfmridir,['c2' segfunc]))
        [c1im, c2im, ~] = segment_funcdat(ppparams.preprocfmridir,[segfunc ',1']);
    else
        c1im = fullfile(ppparams.preprocfmridir,['c1' segfunc]);
        c2im = fullfile(ppparams.preprocfmridir,['c2' segfunc]);
    end

    gmdat = spm_read_vols(spm_vol(c1im));
    wmdat = spm_read_vols(spm_vol(c2im));
    
    braindat = mask;
    braindat((gmdat + wmdat) < 0.2) = 0;
    braindat(braindat > 0.0) = 1;

    % Prepare atlasmaps, componentmaps and brainmask
    [atlasfolder,atlasname,~] = fileparts(params.spatial_networks_masks);
    nettext=fullfile(atlasfolder,[atlasname '.txt']);
    T = readtable(nettext,'FileType','text');

    Vatlas = spm_vol(params.spatial_networks_masks);
    numNetworks = numel(Vatlas);
    atlasVol = spm_read_vols(Vatlas);

    thcompData = reshape(thcompData, [prod(dim(1:3)),numComp]);
    atlasVol = reshape(atlasVol, [prod(dim(1:3)),numNetworks]);
    braindat = reshape(braindat, [prod(dim(1:3)),1]);

    % Spatial sorting
    nvoxbrain = sum(braindat>0);
    useComp = zeros(numNetworks,numComp);

    for inet=1:numNetworks
        regmap = zeros(HInfo.DIM(1:3));
        thregmap = zeros(HInfo.DIM(1:3));

        bpercnetwork(inet) = sum(atlasVol(:,inet)>0)/nvoxbrain;

        for ic = 1:numComp          
            Compdat = thcompData(:,ic);

            totComp = sum(Compdat(Compdat>0),"all");
            brainComp = sum(Compdat(braindat>0),"all");
            brainFract = brainComp/totComp;

            if (brainFract)>0.25
                Compdat(braindat<0.5) = 0;
            
                totComp = sum(Compdat(Compdat>0),"all");
                atlasComp = sum(Compdat(atlasVol(:,inet)>0),"all");
    
                inFract(inet,ic) = atlasComp/totComp;
            
                if inFract(inet,ic)>bpercnetwork(inet)*1.2
                    regmap = regmap+utcompData(:,:,:,ic); 
                    thregmap = thregmap+compData(:,:,:,ic);
    
                    useComp(inet,ic)=1;
                end
            end
        
            clear Compdat
        end

        V = spm_vol(compFiles);

        Vout = V(1);
        Vout.fname = fullfile(ppparams.resultmap,['ica_map_atlas-' atlasname '_network-' num2str(T.Index(inet),'%03d') '-' T.Name{inet} '.nii']);
        Vout = spm_write_vol(Vout,regmap);
    
        Vout = V(1);
        Vout.fname = fullfile(ppparams.resultmap,['thres_ica_map_atlas-' atlasname '_network-' num2str(T.Index(inet),'%03d') '-' T.Name{inet} '.nii']);
        Vout = spm_write_vol(Vout,thregmap);
    
        clear Vout regmap thregmap
    end

    save("components_per_networkMap.mat",'inFract','bpercnetwork','useComp')

end

function [c1im, c2im, c3im] = segment_funcdat(folder,segfunc)

%% Do segmentation of func data

preproc.channel.vols = {fullfile(folder,segfunc)};
preproc.channel.biasreg = 0.001;
preproc.channel.biasfwhm = 60;
preproc.channel.write = [0 0];
preproc.tissue(1).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,1')};
preproc.tissue(1).ngaus = 1;
preproc.tissue(1).native = [1 0];
preproc.tissue(1).warped = [0 0];
preproc.tissue(2).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,2')};
preproc.tissue(2).ngaus = 1;
preproc.tissue(2).native = [1 0];
preproc.tissue(2).warped = [0 0];
preproc.tissue(3).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,3')};
preproc.tissue(3).ngaus = 2;
preproc.tissue(3).native = [1 0];
preproc.tissue(3).warped = [0 0];
preproc.tissue(4).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,4')};
preproc.tissue(4).ngaus = 3;
preproc.tissue(4).native = [0 0];
preproc.tissue(4).warped = [0 0];
preproc.tissue(5).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,5')};
preproc.tissue(5).ngaus = 4;
preproc.tissue(5).native = [0 0];
preproc.tissue(5).warped = [0 0];
preproc.tissue(6).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,6')};
preproc.tissue(6).ngaus = 2;
preproc.tissue(6).native = [0 0];
preproc.tissue(6).warped = [0 0];
preproc.warp.mrf = 1;
preproc.warp.cleanup = 1;
preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
preproc.warp.affreg = 'mni';
preproc.warp.fwhm = 0;
preproc.warp.samp = 3;
preproc.warp.write = [0 1];
preproc.warp.vox = NaN;
preproc.warp.bb = [NaN NaN NaN;NaN NaN NaN];

spm_preproc_run(preproc);

c1im = fullfile(folder,['c1' segfunc]);
c2im = fullfile(folder,['c2' segfunc]);
c3im = fullfile(folder,['c3' segfunc]);
