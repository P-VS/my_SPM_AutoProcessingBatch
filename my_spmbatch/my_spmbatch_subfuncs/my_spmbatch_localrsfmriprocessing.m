function params = my_spmbatch_localrsfmriprocessing(sub,ses,run,task,datpath,params)

%% Search for the data folders

[ppparams,params,datpath] = my_spmbatch_1stlevel_FindData(sub,ses,run,task,datpath,params);

%% Make result map

if params.func.mruns
    ppparams.resultfolder = ['Local-RSfMRI-' task '_' params.analysisname '_run-' num2str(run)];
    ppparams.resultmap = fullfile(ppparams.subpath,['Local-RSfMRI-' task '_' params.analysisname '_run-' num2str(run)]);
else
    ppparams.resultfolder = ['Local-RSfMRI-' task '_' params.analysisname];
    ppparams.resultmap = fullfile(ppparams.subpath,['Local-RSfMRI-' task '_' params.analysisname]);
end

params.resultmap = ppparams.resultfolder;

if ~exist(ppparams.resultmap,'dir'); mkdir(ppparams.resultmap); end

%% fMRI parameters and data

jsondat = fileread(ppparams.frun(1).funcjsonfile);
jsondat = jsondecode(jsondat);

tr = jsondat.RepetitionTime;

Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(1).func(1).funcfile));
nvols = min([numel(Vfunc),50]);
fdata = spm_read_vols(Vfunc);
mask = my_spmbatch_mask(fdata(:,:,:,1:nvols));

Vmask = Vfunc(1);
rmfield(Vmask,'pinfo');
Vmask.fname = fullfile(ppparams.resultmap,['mask_' ppparams.frun(1).func(1).funcfile]);
Vmask.descrip = 'my_spmbatch - mask';
Vmask.dt = [spm_type('float32'),spm_platform('bigend')];
Vmask.n = [1 1];
Vmask = spm_write_vol(Vmask,mask);

ppparams.mask_file = Vmask.fname;

if params.do_ALFF
    fprintf('\nComputing ALFF...');

    fdata = reshape(fdata,[prod(Vfunc(1).dim(1:3)),numel(Vfunc)]);

    alffbrain = my_smbatch_alff(fdata(mask>0,:),params.LF_band,tr);

    alffmap = zeros(Vfunc(1).dim);
    alffmap(mask>0) = alffbrain;
    
    Valff = Vfunc(1);
    rmfield(Valff,'pinfo');
    Valff.fname = fullfile(ppparams.resultmap,['alff_' ppparams.frun(1).func(1).funcfile]);
    Valff.descrip = 'my_spmbatch - alff';
    Valff.dt = [spm_type('float32'),spm_platform('bigend')];
    Valff.n = [1 1];
    Valff = spm_write_vol(Valff,alffmap);

    clear alffmap alffbrain Valff
end

if params.do_fALFF
    fprintf('\nComputing fALFF...');

    fdata = reshape(fdata,[prod(Vfunc(1).dim(1:3)),numel(Vfunc)]);

    alffbrain = my_smbatch_alff(fdata(mask>0,:),params.LF_band,tr);
    ahffbrain = my_smbatch_alff(fdata(mask>0,:),[0,Inf],tr);

    falffmap = zeros(Vfunc(1).dim);
    falffmap(mask>0) = alffbrain./ahffbrain;
    falffmap(isnan(falffmap)) = 0.0;
    
    Vfalff = Vfunc(1);
    rmfield(Vfalff,'pinfo');
    Vfalff.fname = fullfile(ppparams.resultmap,['falff_' ppparams.frun(1).func(1).funcfile]);
    Vfalff.descrip = 'my_spmbatch - falff';
    Vfalff.dt = [spm_type('float32'),spm_platform('bigend')];
    Vfalff.n = [1 1];
    Vfalff = spm_write_vol(Vfalff,falffmap);

    clear alffbrain ahffbrain falffmap Vfalff
end

if params.do_ReHo
    fprintf('\nComputing ReHo...');

    fdata = reshape(fdata,[Vfunc(1).dim(1),Vfunc(1).dim(2),Vfunc(3).dim(1),numel(Vfunc)]);

    rehomap = my_smbatch_reho(fdata,mask,params.Nvoxels);
    
    Vreho = Vfunc(1);
    rmfield(Vreho,'pinfo');
    Vreho.fname = fullfile(ppparams.resultmap,['reho_' ppparams.frun(1).func(1).funcfile]);
    Vreho.descrip = 'my_spmbatch - reho';
    Vreho.dt = [spm_type('float32'),spm_platform('bigend')];
    Vreho.n = [1 1];
    Vreho = spm_write_vol(Vreho,rehomap);

    clear rehomap Vreho
end

clear fdata Vfunc