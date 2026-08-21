function my_spmbatch_conn_GroupSBC(sublist,nsessions,datpath,params)

if ~params.func.mruns, params.func.runs = [1]; end

ppparams.resultfolder = fullfile(params.outfolder,params.analysisname);

if ~exist(ppparams.resultfolder,'dir'); mkdir(ppparams.resultfolder); end

clear batch;
batch.filename=fullfile(ppparams.resultfolder,['conn_' params.analysisname '.mat']);

%% BATCH.parallel DEFINES PARALLELIZATION OPTIONS (applies to any Setup/Setup.preprocessing/Denoising/Analysis/QA steps) %!

if params.use_parallel
    batch.parallel.N = params.maxprocesses;
    batch.parallel.immediatereturn = 0;
else 
    batch.parallel.N = 0;
end

%% BATCH.Setup DEFINES EXPERIMENT SETUP AND PERFORMS INITIAL DATA EXTRACTION AND/OR PREPROCESSING STEPS %!

if exist(batch.filename,'file') && params.add
    batch.Setup.isnew = 0;
    batch.Setup.overwrite = 0;
    batch.Setup.add = 1;
else
    batch.Setup.isnew = 1;
    batch.Setup.overwrite = 1;
end

batch.Setup.done = 1; 

batch.Setup.nsubjects = numel(sublist);

ilist = 0;
filelist = {};

nses = numel(nsessions)*numel(params.task)*numel(params.func.runs);

for ic=1:nses
    for ip=1:batch.Setup.nsubjects
        for is=1:nses
            batch.Setup.conditions.onsets{ic}{ip}{is}=[];
            batch.Setup.conditions.durations{ic}{ip}{is}=[];
        end
    end
end

for isub=1:batch.Setup.nsubjects
    for isess=1:numel(nsessions)
        for itask=1:numel(params.task)
            for irun=1:numel(params.func.runs)
                substring = ['sub-' num2str(sublist(isub),['%0' num2str(params.sub_digits) 'd'])];
                sesstring = ['ses-' num2str(nsessions(isess),'%03d')];
        
                subfolder = fullfile(datpath,substring,sesstring);
        
                % Search functional data
                preprocfmridir = fullfile(subfolder,params.preprocfmridir);

                funcfile = myspmb_connSBC_findDat(preprocfmridir,substring,sesstring,params.task{itask},params.func.runs(irun),params,'func');
        
                curses = ((isess-1)*numel(params.task)+(itask-1))*numel(params.func.runs)+irun;

                if ~isempty(funcfile), batch.Setup.functionals{isub}{curses} = funcfile; end

                ilist = ilist+1;
                filelist{ilist}.conn_subid = ['sub-' num2str(isub,['%0' num2str(params.sub_digits) 'd'])];
                filelist{ilist}.conn_sessid = ['sess-' num2str(curses,'%02d')];
                filelist{ilist}.conn_scanid = [substring '/' sesstring '/task-' params.task{itask} '/run-' num2str(params.func.runs(irun))];
                filelist{ilist}.file = funcfile;
            
                % Search TR
                if curses==1
                    funcdir = fullfile(subfolder,'func');
                
                    jsonfile = myspmb_connSBC_findDat(funcdir,substring,sesstring,params.task{1},params.func.runs(1),params,'json');
                
                    jsondat = fileread(jsonfile);
                    jsondat = jsondecode(jsondat);
                
                    tr = jsondat.RepetitionTime;
                
                    batch.Setup.RT(isub) = tr;
                end

                % Search task onsets and durations
                tsvfile = myspmb_connSBC_findDat(funcdir,substring,sesstring,params.task{1},params.func.runs(1),params,'tsv');
            
                if ~isempty(tsvfile)
                    batch.Setup.conditions.importfile{isub}{curses} = tsvfile;
                else
                    condname = ['rest_task-' params.task{itask} '_ses-' num2str(nsessions(isess),'%03d') '_run-' num2str(params.func.runs(irun))];
                    if isfield(batch.Setup.conditions,'names')
                        if ~any(strcmp(condname,batch.Setup.conditions.names))
                            batch.Setup.conditions.names = {batch.Setup.conditions.names{:} condname};
                        end
                    else
                        batch.Setup.conditions.names = {condname};
                    end
                    tmp = find(strcmp(condname,batch.Setup.conditions.names));
                    batch.Setup.conditions.onsets{tmp(1)}{isub}{curses}=0;
                    batch.Setup.conditions.durations{tmp(1)}{isub}{curses}=Inf;
                end

                % Search anatomical data
                anatdir = fullfile(subfolder,'preproc_anat');
            
                anatfile = myspmb_connSBC_findDat(anatdir,substring,sesstring,'','',params,'anat');
                batch.Setup.structurals{isub}{curses} = anatfile;
            
                gmfile = myspmb_connSBC_findDat(anatdir,substring,sesstring,'','',params,'GM');
                if ~isempty(gmfile), batch.Setup.masks.Grey{isub}{curses} = gmfile; else batch.Setup.preprocessing.steps= {'structural_segment'}; end
            
                wmfile = myspmb_connSBC_findDat(anatdir,substring,sesstring,'','',params,'WM');
                if ~isempty(wmfile), batch.Setup.masks.White{isub}{curses} = wmfile; else batch.Setup.preprocessing.steps= {'structural_segment'}; end
                
                csffile = myspmb_connSBC_findDat(anatdir,substring,sesstring,'','',params,'CSF');
                if ~isempty(csffile), batch.Setup.masks.CSF{isub}{curses} = csffile; else batch.Setup.preprocessing.steps= {'structural_segment'}; end
            end
        end
    end
end
save('setup.mat','batch')
Setup.secondarydatasets.functionals_type = 2; %same files as functional data field after removing leading 's' from filename

batch.Setup.rois.names={'atlas'};
batch.Setup.rois.files{1}=fullfile(fileparts(which('conn')),'rois','atlas.nii');
batch.Setup.rois.dimensions = [1];
batch.Setup.rois.multiplelabels = [1];

for ic=1:numel(params.sub_covariates)
    batch.Setup.subjects.effects{ic} = params.sub_covariates{ic}.vector;
    batch.Setup.subjects.effect_names{ic} = params.sub_covariates{ic}.names;
end

batch.Setup.analyses = [];
if params.do_ROI2ROI, batch.Setup.analyses = [batch.Setup.analyses 1]; end
if params.do_Seed2voxel, batch.Setup.analyses = [batch.Setup.analyses 2]; end
if params.do_voxel2voxel, batch.Setup.analyses = [batch.Setup.analyses 3]; end

batch.Setup.outputfiles = [1,0,1,1,0,0]; %Optional output files (outputfiles(1): 1/0 creates confound beta-maps; outputfiles(2): 1/0 creates 
%                                         confound-corrected timeseries; outputfiles(3): 1/0 creates seed-to-voxel r-maps) ;outputfiles(4): 
%                                         1/0 creates seed-to-voxel p-maps) ;outputfiles(5): 1/0 creates seed-to-voxel FDR-p-maps); 
%                                         outputfiles(6): 1/0 creates ROI-extraction REX files;

%% CONN Denoising
batch.Denoising.done=1;
batch.Denoising.filter=[0.01, 0.1];          % frequency filter (band-pass values, in Hz)
batch.Denoising.detrending = 1;
batch.Denoising.despiking = 0;
batch.Denoising.regbp = 1;
for ic=1:numel(batch.Setup.conditions.names)
    batch.Denoising.confounds{ic} = ['Effect of ' batch.Setup.conditions.names{ic}];
end

%% CONN Analysis
% BATCH.Analysis PERFORMS FIRST-LEVEL ANALYSES (ROI-to-ROI and seed-to-voxel) %!
if params.do_ROI2ROI || params.do_Seed2voxel
    batch.Analysis.done=1;
    batch.Analysis.overwrite = 1;
    batch.Analysis.name = 'Correlation All ROIs'; 
    if params.do_ROI2ROI, batch.Analysis.type = 1; else batch.Analysis.type = 0; end
    if params.do_Seed2voxel, batch.Analysis.type = batch.Analysis.type + 2; end
    batch.Analysis.measure = 1;               % connectivity measure used {1 = 'correlation (bivariate)', 2 = 'correlation (semipartial)', 3 = 'regression (bivariate)', 4 = 'regression (multivariate)';
    batch.Analysis.weight = 2;                % within-condition weight used {1 = 'none', 2 = 'hrf', 3 = 'hanning';
    batch.Analysis.sources = {'atlas'};              % (defaults to all ROIs)
end

conn_batch(batch);

%% -------------------------------------
function file = myspmb_connSBC_findDat(subfolder,substring,sesstring,task,run,params,dattype)

file = '';

if ~isfolder(subfolder), return; end

namefilters(1).name = substring;
namefilters(1).required = true;

namefilters(2).name = sesstring;
namefilters(2).required = false;

if contains(dattype,'func') || contains(dattype,'json')
    namefilters(3).name = ['task-' task];
    namefilters(3).required = true;

    namefilters(4).name = ['run-' num2str(run)];
    namefilters(4).required = params.func.mruns;

    namefilters(5).name = '_bold';
    namefilters(5).required = true;
    
    if contains(dattype,'func')
        namefilters(6).name = params.fmri_prefix;
        namefilters(6).required = true;
    end
elseif contains(dattype,'anat')
    namefilters(3).name = ['_T1w'];
    namefilters(3).required = true;
elseif contains(dattype,'tsv')
    namefilters(3).name = ['task-' task];
    namefilters(3).required = true;

    namefilters(4).name = ['run-' num2str(run)];
    namefilters(4).required = params.func.mruns;

    namefilters(5).name = '_events';
    namefilters(5).required = true;
end

if contains(dattype,'json') 
    jsonlist = my_spmbatch_dirfilelist(subfolder,'json',namefilters,false);

    if isempty(jsonlist), return; end
elseif contains(dattype,'tsv') 
    tsvlist = my_spmbatch_dirfilelist(subfolder,'tsv',namefilters,false);

    if isempty(tsvlist), return; end
else
    niilist = my_spmbatch_dirfilelist(subfolder,'nii',namefilters,false);

    if isempty(niilist), return; end
end

switch dattype
    case 'func'
        prefixlist = split({niilist.name},'sub-');
        if numel(niilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end
        
        tmp = find(strcmp(prefixlist,params.fmri_prefix));
        file = fullfile(subfolder,niilist(tmp(1)).name);
    case 'json'
        file = fullfile(subfolder,jsonlist(1).name);
    case 'anat'
        tmp = find(contains({niilist.name},'_Crop_1'));
        if ~isempty(tmp), niilist = niilist(tmp); end

        prefixlist = split({niilist.name},'sub-');
        prefixlist = prefixlist(:,:,1);

        tmp = find(strcmp(prefixlist,'we'));
        if isempty(tmp), tmp = find(strcmp(prefixlist,'wme')); end

        if ~isempty(tmp), file = fullfile(subfolder,niilist(tmp).name); end
    case 'GM'
        tmp = find(contains({niilist.name},'_Crop_1'));
        if ~isempty(tmp), niilist = niilist(tmp); end

        prefixlist = split({niilist.name},'sub-');
        prefixlist = prefixlist(:,:,1);

        tmp = find(strcmp(prefixlist,'wc1e'));
        if isempty(tmp), tmp = find(strcmp(prefixlist,'wmp1e')); end

        if ~isempty(tmp), file = fullfile(subfolder,niilist(tmp).name); end
    case 'WM'
        tmp = find(contains({niilist.name},'_Crop_1'));
        if ~isempty(tmp), niilist = niilist(tmp); end

        prefixlist = split({niilist.name},'sub-');
        prefixlist = prefixlist(:,:,1);

        tmp = find(strcmp(prefixlist,'wc2e'));
        if isempty(tmp), tmp = find(strcmp(prefixlist,'wmp2e')); end

        if ~isempty(tmp), file = fullfile(subfolder,niilist(tmp).name); end
    case 'CSF'
        tmp = find(contains({niilist.name},'_Crop_1'));
        if ~isempty(tmp), niilist = niilist(tmp); end

        prefixlist = split({niilist.name},'sub-');
        prefixlist = prefixlist(:,:,1);

        tmp = find(strcmp(prefixlist,'wc3e'));
        if isempty(tmp), tmp = find(strcmp(prefixlist,'wmp3e')); end

        if ~isempty(tmp), file = fullfile(subfolder,niilist(tmp).name); end
    case 'tsv'
        file = fullfile(subfolder,tsvlist(1).name);
end


%% -------------------------------------
function spmmat_file = myspmb_connSBC_makeSPMmat(subfolder,substring,sesstring,task,run,params,dattype)
