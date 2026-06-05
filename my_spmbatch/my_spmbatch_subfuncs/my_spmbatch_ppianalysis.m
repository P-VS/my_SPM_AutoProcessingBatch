function params = my_spmbatch_ppianalysis(sub,ses,run,task,datpath,params)

%% Search for the data folders

ppparams.substring = ['sub-' num2str(sub)];
if isfolder(fullfile(datpath,['sub-' num2str(sub)])), ppparams.substring = ['sub-' num2str(sub)]; end
if isfolder(fullfile(datpath,['sub-0' num2str(sub)])), ppparams.substring = ['sub-0' num2str(sub)]; end
if isfolder(fullfile(datpath,['sub-00' num2str(sub)])), ppparams.substring = ['sub-00' num2str(sub)]; end
if isfolder(fullfile(datpath,['sub-000' num2str(sub)])), ppparams.substring = ['sub-000' num2str(sub)]; end

ppparams.sesstring = ['ses-' num2str(ses)];
if isfolder(fullfile(datpath,ppparams.substring,['ses-' num2str(ses)])), ppparams.sesstring = ['ses-' num2str(ses)]; end
if isfolder(fullfile(datpath,ppparams.substring,['ses-0' num2str(ses)])), ppparams.sesstring = ['ses-0' num2str(ses)]; end
if isfolder(fullfile(datpath,ppparams.substring,['ses-00' num2str(ses)])), ppparams.sesstring = ['ses-00' num2str(ses)]; end

ppparams.subpath = fullfile(datpath,ppparams.substring,ppparams.sesstring);

if ~isfolder(ppparams.subpath), ppparams.subpath = fullfile(datpath,ppparams.substring); end

if ~isfolder(ppparams.subpath)
    fprintf(['No data folder for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.fmriresultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.SPMMAT_analysisname]);

if ~isfolder(ppparams.fmriresultmap)
    fprintf(['No SPMMAT folder for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.ppiresultmap = fullfile(ppparams.fmriresultmap,['PPI-' task '_' params.PPI_analysisname]);
if exist(ppparams.ppiresultmap,'dir'); rmdir(ppparams.ppiresultmap,'s'); end
mkdir(ppparams.ppiresultmap)

tmfc.project_path = ppparams.ppiresultmap;

tmfc.defaults.parallel = 0;      
tmfc.defaults.maxmem = 2^32;
tmfc.defaults.resmem = true;
tmfc.defaults.analysis = 1;

%% Import SPM.mat data

load(fullfile(ppparams.fmriresultmap,'SPM.mat'));

% SPM output directory
SPM.swd = ppparams.fmriresultmap;

P = SPM.xY.P;
rmfield(SPM.xY,'P');

% SPM fMRI data files
for ifile=1:size(P,1)
    [ffolder,fname,fext] = fileparts(P(ifile,:));

    if ifile==1, SPM.xY.P = fullfile(ppparams.subpath,params.preprocfmridir,[fname fext]); else SPM.xY.P(ifile,:) = fullfile(ppparams.subpath,params.preprocfmridir,[fname fext]); end
    SPM.xY.VY(ifile).fname = fullfile(ppparams.subpath,params.preprocfmridir,[fname '.nii']);
end

[ffolder,fname,fext] = fileparts(SPM.xM.VM.fname);
SPM.xM.VM.fname = fullfile(ppparams.subpath,params.preprocfmridir,[fname fext]);

save(fullfile(ppparams.ppiresultmap,'SPM.mat'),'SPM')

spmfiles = dir(ppparams.fmriresultmap);
confiles = find(startsWith({spmfiles.name},'con_'));
for iconf=1:numel(confiles)
    copyfile(fullfile(ppparams.fmriresultmap,spmfiles(confiles(iconf)).name),fullfile(ppparams.ppiresultmap,spmfiles(confiles(iconf)).name));
end
betafiles = find(startsWith({spmfiles.name},'beta_'));
for iconf=1:numel(betafiles)
    copyfile(fullfile(ppparams.fmriresultmap,spmfiles(betafiles(iconf)).name),fullfile(ppparams.ppiresultmap,spmfiles(betafiles(iconf)).name));
end
copyfile(fullfile(ppparams.fmriresultmap,'mask.nii'),fullfile(ppparams.ppiresultmap,'mask.nii'));
copyfile(fullfile(ppparams.fmriresultmap,'RPV.nii'),fullfile(ppparams.ppiresultmap,'RPV.nii'));
copyfile(fullfile(ppparams.fmriresultmap,'ResMS.nii'),fullfile(ppparams.ppiresultmap,'ResMS.nii'));

tmfc.subjects(1).path = fullfile(ppparams.ppiresultmap,'SPM.mat');
tmfc.subjects(1).name = ppparams.substring;

%% Import ROI data
tmfc.ROI_set_number = 1;
tmfc.ROI_set(1).set_name = params.PPI_analysisname;
tmfc.ROI_set(1).type = 'binary_images';

Vmask = spm_vol(fullfile(SPM.swd,'mask.nii'));
mask = spm_read_vols(Vmask);
roi_path = fullfile(ppparams.ppiresultmap,'ROI_sets',params.PPI_analysisname,'masked_ROIs');

if ~isfolder(roi_path), mkdir(roi_path); end

roilist = dir(params.VOIfolder);
if isempty(roilist), return; end

tmp = find(strlength({roilist.name})>4); %Remove '.' and '..'
if ~isempty(tmp)
    roilist = roilist(tmp);
else
    return 
end

tmp = find(~contains({roilist.name},'._')); %Remove the hiden files from Mac from the list
if ~isempty(tmp)
    roilist = roilist(tmp);
else
    return 
end

tmp = find(contains({roilist.name},'.nii'));%Filter list based on the file extension
if ~isempty(tmp) 
    roilist = roilist(tmp); 
else 
    return 
end

for iroi=1:numel(roilist)
    ROI = spm_vol(fullfile(roilist(iroi).folder,roilist(iroi).name));

    if ~(ROI.dim==Vmask.dim)
        fprintf(['ROI ' roilist(iroi).name ' not of the same dimensions as the subject mask'])
        fprintf('\nPP_Error\n');
        return
    end

    ROImask = spm_read_vols(ROI);

    ROImask = (ROImask.*mask)>0.1;

    ROI.fname = fullfile(roi_path,roilist(iroi).name);
    ROI = spm_write_vol(ROI,ROImask);

    Roiname = split(roilist(iroi).name,'.nii');
    tmfc.ROI_set(1).ROIs(iroi).name = Roiname{1};
    tmfc.ROI_set(1).ROIs(iroi).path_masked =ROI.fname;

    if numel(ROImask>0.1)==0
        fprintf(['ROI ' roilist(iroi).name ' no voxels in ROI'])
        fprintf('\nPP_Error\n');
        return
    end
end

%% VOI
for icond=1:numel(SPM.Sess(1).U)
    tmfc.ROI_set(1).gPPI.conditions(icond).sess   = 1;   
    tmfc.ROI_set(1).gPPI.conditions(icond).number = icond; 
    tmfc.ROI_set(1).gPPI.conditions(icond).pmod   = 1; 
    tmfc.ROI_set(1).gPPI.conditions(icond).name = SPM.Sess(1).U(icond).name; 
    tmfc.ROI_set(1).gPPI.conditions(icond).file_name = ['[Sess_1]_[Cond_' num2str(icond) ']_[' SPM.Sess(1).U(icond).name{1} ']'];
end

disp('Initiating VOI computation...');
sub_check = tmfc_VOI(tmfc,1,1); 
tmfc.ROI_set(1).subjects(1).VOI = sub_check(1);
disp('VOI computation completed.');

save(fullfile(tmfc.project_path,'tmfc_autosave.mat'),'tmfc');

if params.doPPI
    %% PPI
    disp('Initiating PPI computation...');
             
    % Define mean centering
    centering = 'with_mean_centering';
    
    tmfc.ROI_set(1).PPI_centering = centering;
    
    % Run PPI calculation
    sub_check = tmfc_PPI(tmfc,tmfc.ROI_set_number,1);
    tmfc.ROI_set(1).subjects(1).PPI = sub_check(1);
    
    disp('PPI computation completed.');
    
    save(fullfile(tmfc.project_path,'tmfc_autosave.mat'),'tmfc');
    
    %% gPPI  
    disp('Initiating gPPI computation...');
    
    [sub_check, contrasts] = tmfc_gPPI(tmfc,1,1);    
    tmfc.ROI_set(1).subjects(1).gPPI = sub_check(1);
    
    for iCon = 1:length(contrasts)
        tmfc.ROI_set(1).contrasts.gPPI(iCon).title = contrasts(iCon).title;
        tmfc.ROI_set(1).contrasts.gPPI(iCon).weights = contrasts(iCon).weights;
    end
    
    disp('gPPI computation completed.');

    save(fullfile(tmfc.project_path,'tmfc_autosave.mat'),'tmfc');
end

if params.doBSC
    %% LSS GLM
    disp('Initiating LSS computation...');

    tmfc.LSS.conditions = tmfc.ROI_set(1).gPPI.conditions;

    sub_check = tmfc_LSS(tmfc,1);  
    tmfc.subjects(1).LSS = sub_check(1);

    save(fullfile(tmfc.project_path,'tmfc_autosave.mat'),'tmfc');
    disp('LSS computation completed.');

    %% BSC
    disp('Initiating BSC LSS computation...');  

    tmfc.ROI_set(1).BSC = 'first_eigenvariate';

    [sub_check, contrasts] = tmfc_BSC(tmfc,1);
    tmfc.ROI_set(1).subjects(1).BSC = sub_check(1);

    for iCon = 1:length(contrasts)
        tmfc.ROI_set(1).contrasts.BSC(iCon).title = contrasts(iCon).title;
        tmfc.ROI_set(1).contrasts.BSC(iCon).weights = contrasts(iCon).weights;
    end
    save(fullfile(tmfc.project_path,'tmfc_autosave.mat'),'tmfc');
    disp('BSC LSS computation completed.');
end