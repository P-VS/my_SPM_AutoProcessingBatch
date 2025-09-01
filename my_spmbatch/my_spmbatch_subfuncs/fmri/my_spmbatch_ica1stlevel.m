function params = my_spmbatch_ica1stlevel(sub,ses,run,task,datpath,params)

curdir = pwd;

matlabbatch = {};

%% Search for the data folders
ppparams = my_spmbatch_1stlevel_FindData(sub,ses,run,task,datpath,params);

%% fMRI model specification
if params.func.mruns && contains(params.func.use_runs,'separately')
    ppparams.resultmap = fullfile(ppparams.subpath,['ICA-' task '_' params.analysisname '_run-' num2str(run)]);
else
    ppparams.resultmap = fullfile(ppparams.subpath,['ICA-' task '_' params.analysisname]);
end

if exist(ppparams.resultmap,'dir'); rmdir(ppparams.resultmap,'s'); end
mkdir(ppparams.resultmap)

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

%% Adapt template parameter to specific subject data
% Get t_r
jsondat = jsondecode(fileread(ppparams.frun(1).funcjsonfile));
t_r = jsondat.("RepetitionTime");

for ir=1:numel(params.iruns)
    fsplit = split(ppparams.ppfmridat{ir}.sess{1}.func{1,1},'.nii');
    ica_source_file{ir} = [fsplit{1} '.nii'];
end

icatb_defaults;

% Load template parameters .mat file
icaparms_file = ica_1stlevel_make_gift_parameters(ica_source_file,t_r,numel(params.iruns),params,ppparams);

%% Set up ICA
icaparam_file = icatb_setup_analysis(icaparms_file);

load(icaparam_file);

%% Run Analysis (All steps)
sesInfo = icatb_runAnalysis(sesInfo, [2:6]);

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

compData = spm_read_vols(spm_vol(compFiles));
dim = size(compData);
HInfo.V = spm_vol(structFile);
HInfo.DIM = dim(1:3);
HInfo.VOX = double(HInfo.V(1).private.hdr.pixdim(2:4)); HInfo.VOX = abs(HInfo.VOX);

% load time course
icaTimecourse = icatb_loadICATimeCourse(compFiles, 'real', [], [1:numComp]);

% Structural volume
dim = HInfo.DIM;
tdim = size(icaTimecourse(:,1));

%% Load model timecourses from designa matrix
load(SPM_file)

modelTimecourse = SPM.xX.X;
szX = size(modelTimecourse);

refInfo.spmMatFlag = 1;
refInfo.SPMFile = {SPM_file};
for ireg=1:szX(2)
    refInfo.selectedRegressors(ireg).name = SPM.xX.name{ireg};
end
refInfo.modelIndex = [1:szX(2)];

%% Temporal sorting components

[sortParameters] = icatb_sortComponents('sortingCriteria', 'multiple regression', ...
        'sortingType', 'temporal', ...
        'icaTimecourse', icaTimecourse, ...
        'modelTimecourse', modelTimecourse, ...
        'num_Regress', szX(2), ...
        'num_DataSets', 1, ...
        'refInfo', refInfo, ...
        'diffTimePoints', szX(1), ...
        'numcomp', numComp, ...
        'icasig', [], ...
        'structHInfo', HInfo, ...
        'num_sort_subjects',  1, ...
        'num_sort_sessions', 1, ...
        'viewingSet', 'All datasets', ...
        'input_prefix', 'ica_', ...
        'output_dir', ppparams.resultmap);

T = readtable(fullfile(ppparams.resultmap,'ica__temporal_partial_corr.txt'),NumHeaderLines=2,Delimiter="\t");

sComps = table2array(T(1,2:end));
regNames = table2array(T(2:(szX(2)+1),1));
partCorrComps = table2array(T(2:(szX(2)+1),2:end));

V = spm_vol(compFiles);
for ireg=1:szX(2)-1
    regmap = zeros(HInfo.DIM(1:3));
    for ic=1:numComp
        if partCorrComps(ireg,ic)>0.4, regmap = regmap+compData(:,:,:,sComps(ic)); end
    end

    thregmap = reshape(regmap, 1, prod(HInfo.DIM(1:3)));
    [thregmap] = icatb_applyDispParameters(thregmap, 1, 3, 1, dim(1:3), HInfo);
    thregmap = reshape(thregmap, [dim(1), dim(2), dim(3)]);

    sregName = split(regNames{ireg},'Sn(1) ');

    Vout = V(1);
    Vout.fname = fullfile(ppparams.resultmap,['ica_map_' sregName{end} '.nii']);
    Vout = spm_write_vol(Vout,regmap);

    Vout = V(1);
    Vout.fname = fullfile(ppparams.resultmap,['thres_ica_map_' sregName{end} '.nii']);
    Vout = spm_write_vol(Vout,thregmap);

    clear Vout regmap thregmap
end

%% Contrasts

for ic=1:numel(params.contrast)
    contDat = zeros(HInfo.DIM(1:3));
    conName = 'ica_map';

    for icond=1:numel(params.contrast(ic).conditions)
        tmp=find(contains(regNames,params.contrast(ic).conditions{icond}));
        if ~isempty(tmp)
            sregName = split(regNames{tmp(1)},'Sn(1) ');
            regMap = fullfile(ppparams.resultmap,['ica_map_' sregName{end} '.nii']);

            regDat = spm_read_vols(spm_vol(regMap));

            contDat = contDat+params.contrast(ic).vector(icond)*regDat;

            clear regDat

            if params.contrast(ic).vector(icond)>0; conName = [conName '_Pos-' params.contrast(ic).conditions{icond}]; end
            if params.contrast(ic).vector(icond)<0; conName = [conName '_Neg-' params.contrast(ic).conditions{icond}]; end
        end
    end

    V = spm_vol(compFiles);

    Vout = V(1);
    Vout.fname = fullfile(ppparams.resultmap,[conName '.nii']);
    Vout = spm_write_vol(Vout,contDat);

    clear contDat Vout
end

cd(curdir)