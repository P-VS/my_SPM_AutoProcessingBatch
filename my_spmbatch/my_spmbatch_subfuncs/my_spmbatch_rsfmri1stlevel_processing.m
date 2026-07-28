function params = my_spmbatch_rsfmri1stlevel_processing(sub,ses,run,task,datpath,params)

%% Search for the data folders

[ppparams,params,datpath] = my_spmbatch_1stlevel_FindData(sub,ses,run,task,datpath,params);

ppparams.preproc_anat = fullfile(ppparams.subpath,'preproc_anat');

namefilters(1).name = ppparams.substring;
namefilters(1).required = true;

namefilters(2).name = ppparams.sesstring;
namefilters(2).required = false;

namefilters(3).name = ['_T1w'];
namefilters(3).required = true;

ppanatniilist = my_spmbatch_dirfilelist(ppparams.preproc_anat,'nii',namefilters,false);

if ~isempty(ppanatniilist)
    tmp = find(contains({ppanatniilist.name},'_Crop_1'));
    if ~isempty(tmp), ppanatniilist = ppanatniilist(tmp); end

    prefixlist = split({ppanatniilist.name},'sub-');
    prefixlist = prefixlist(:,:,1);
end

if ~isempty(ppanatniilist)
    tmp = find(strcmp(prefixlist,'we'));

    if isempty(tmp), tmp = find(strcmp(prefixlist,'wme')); end

    if ~isempty(tmp), ppparams.wsubanat = ppanatniilist(tmp).name; end
else
    ppparams.wsubanat = '';
end

%% Make result map

if params.func.mruns
    ppparams.resultfolder = ['CONN-RSfMRI-' task '_' params.analysisname '_run-' num2str(run)];
    ppparams.connmap = fullfile(ppparams.subpath,['CONN-RSfMRI-' task '_' params.analysisname '_run-' num2str(run)]);
    ppparams.resultmap = fullfile(ppparams.connmap,['SPMMAT-RSfMRI-' task '_' params.analysisname '_run-' num2str(run)]);
else
    ppparams.resultfolder = ['CONN-RSfMRI-' task '_' params.analysisname];
    ppparams.connmap = fullfile(ppparams.subpath,['CONN-RSfMRI-' task '_' params.analysisname]);
    ppparams.resultmap = fullfile(ppparams.connmap,['SPMMAT-RSfMRI-' task '_' params.analysisname]);
end

params.resultmap = ppparams.resultfolder;

if ~exist(ppparams.connmap,'dir'); mkdir(ppparams.connmap); end
if exist(ppparams.resultmap,'dir'); rmdir(ppparams.resultmap,'s'); end
mkdir(ppparams.resultmap)

[fmri_spec,ppparams] = my_spmbatch_1stlevel_DefineRSModel(sub,ses,run,task,datpath,params,ppparams);

spm_run_fmri_spec(fmri_spec);

jsondat = fileread(ppparams.frun(1).funcjsonfile);
jsondat = jsondecode(jsondat);

tr = jsondat.RepetitionTime;

batch.filename=fullfile(ppparams.connmap,['conn_singlesubject_nopreproc_' ppparams.substring '_' task '-' params.analysisname '.mat']);
batch.Setup.spmfiles=cellstr(fullfile(ppparams.resultmap,'SPM.mat'));
batch.Setup.structurals=cellstr(fullfile(ppparams.preproc_anat,ppparams.wsubanat));
batch.Setup.nsubjects=1;
batch.Setup.RT=tr;
batch.Setup.rois.names={'atlas'};
batch.Setup.rois.files{1}=fullfile(fileparts(which('conn')),'rois','atlas.nii');
batch.Setup.isnew=1;
batch.Setup.done=1;

%% CONN Denoising
batch.Denoising.filter=[0.01, 0.1];          % frequency filter (band-pass values, in Hz)
batch.Denoising.done=1;

%% CONN Analysis
batch.Analysis.analysis_number=1;       % Sequential number identifying each set of independent first-level analyses
batch.Analysis.measure=1;               % connectivity measure used {1 = 'correlation (bivariate)', 2 = 'correlation (semipartial)', 3 = 'regression (bivariate)', 4 = 'regression (multivariate)';
batch.Analysis.weight=2;                % within-condition weight used {1 = 'none', 2 = 'hrf', 3 = 'hanning';
batch.Analysis.sources={};              % (defaults to all ROIs)
batch.Analysis.done=1;

conn_batch(batch);