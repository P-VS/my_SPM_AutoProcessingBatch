function my_spmbatch_conn_GroupSBC(sublist,nsessions,datpath,params)

if ~params.func.mruns, params.func.runs = [1]; end

ppparams.resultfolder = fullfile(params.outfolder,params.analysisname);

if ~exist(ppparams.resultfolder,'dir'); mkdir(ppparams.resultfolder); end

clear batch;
batch.filename=fullfile(ppparams.resultfolder,['conn_' params.analysisname '.mat']);

%% BATCH.parallel DEFINES PARALLELIZATION OPTIONS (applies to any Setup/Setup.preprocessing/Denoising/Analysis/QA steps) %!

if params.use_parallel
    batch.parallel.N = params.maxprocesses;
    if params.run_background, batch.parallel.immediatereturn = 1; else batch.parallel.immediatereturn = 0; end
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

batch.Setup.nsubjects = numel(sublist);

ilist = 0;
filelist = {};

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
            end
        end
    end
    sesstring = 'ses-001';

    % Search TR
    funcdir = fullfile(subfolder,'func');

    jsonfile = myspmb_connSBC_findDat(funcdir,substring,sesstring,params.task{1},params.func.runs(1),params,'json');

    jsondat = fileread(jsonfile);
    jsondat = jsondecode(jsondat);

    tr = jsondat.RepetitionTime;

    batch. Setup.RT(isub) = tr;

    % Search anatomical data
    anatdir = fullfile(subfolder,'preproc_anat');

    anatfile = myspmb_connSBC_findDat(anatdir,substring,sesstring,'','',params,'anat');
    batch.Setup.structurals{isub} = anatfile;

    gmfile = myspmb_connSBC_findDat(anatdir,substring,sesstring,'','',params,'GM');
    batch.Setup.masks.Grey{isub} = gmfile;

    wmfile = myspmb_connSBC_findDat(anatdir,substring,sesstring,'','',params,'WM');
    batch.Setup.masks.White{isub} = wmfile;
    
    csffile = myspmb_connSBC_findDat(anatdir,substring,sesstring,'','',params,'CSF');
    batch.Setup.masks.CSF{isub} = csffile;
end

Setup.secondarydatasets.functionals_type = 2; %same files as functional data field after removing leading 's' from filename

batch.Setup.rois.names={'atlas'};
batch.Setup.rois.files{1}=fullfile(fileparts(which('conn')),'rois','atlas.nii');


%%-------------------------------------
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
    
    namefilters(6).name = params.fmri_prefix;
    namefilters(6).required = true;
else
    namefilters(3).name = ['_T1w'];
    namefilters(3).required = true;
end

if ~contains(dattype,'json')
    niilist = my_spmbatch_dirfilelist(subfolder,'nii',namefilters,false);

    if isempty(niilist), return; end
else 
    jsonlist = my_spmbatch_dirfilelist(subfolder,'json',namefilters,false);

    if isempty(jsonlist), return; end
end

switch dattype
    case 'func'
        prefixlist = split({niilist.name},'sub-');
        if numel(niilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end
        
        tmp = find(strcmp(prefixlist,params.fmri_prefix));
        file = fullfile(subfolder,niilist(tmp(1)).name);
    case 'json'
        file = fullfile(subfolder,funcjsonlist(1).name);
    case 'anat'
        tmp = find(contains({niilist.name},'_Crop_1'));
        if ~isempty(tmp), niilist = niilist(tmp); end

        prefixlist = split({niilist.name},'sub-');
        prefixlist = prefixlist(:,:,1);

        tmp = find(strcmp(prefixlist,'we'));
        if isempty(tmp), tmp = find(strcmp(prefixlist,'wme')); end

        if ~isempty(tmp), file = niilist(tmp).name; end
    case 'GM'
        tmp = find(contains({niilist.name},'_Crop_1'));
        if ~isempty(tmp), niilist = niilist(tmp); end

        prefixlist = split({niilist.name},'sub-');
        prefixlist = prefixlist(:,:,1);

        tmp = find(strcmp(prefixlist,'wc1e'));
        if isempty(tmp), tmp = find(strcmp(prefixlist,'wmp1e')); end

        if ~isempty(tmp), file = niilist(tmp).name; end
    case 'WM'
        tmp = find(contains({niilist.name},'_Crop_1'));
        if ~isempty(tmp), niilist = niilist(tmp); end

        prefixlist = split({niilist.name},'sub-');
        prefixlist = prefixlist(:,:,1);

        tmp = find(strcmp(prefixlist,'wc2e'));
        if isempty(tmp), tmp = find(strcmp(prefixlist,'wmp2e')); end

        if ~isempty(tmp), file = niilist(tmp).name; end
    case 'CSF'
        tmp = find(contains({niilist.name},'_Crop_1'));
        if ~isempty(tmp), niilist = niilist(tmp); end

        prefixlist = split({niilist.name},'sub-');
        prefixlist = prefixlist(:,:,1);

        tmp = find(strcmp(prefixlist,'wc3e'));
        if isempty(tmp), tmp = find(strcmp(prefixlist,'wmp3e')); end

        if ~isempty(tmp), file = niilist(tmp).name; end
end