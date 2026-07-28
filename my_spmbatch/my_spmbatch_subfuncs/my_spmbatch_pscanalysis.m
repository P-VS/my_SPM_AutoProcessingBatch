function params = my_spmbatch_pscanalysis(sub,ses,run,task,datpath,params)

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
    e.message = ['No data folder for subject ' num2str(sub) ' session ' num2str(ses)];
    error(e)
    return
end

if params.func.mruns
    ppparams.resultfolder = ['SPMMAT-' task '_' params.SPMMAT_analysisname '_run-' num2str(run)];
    ppparams.fmriresultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.SPMMAT_analysisname '_run-' num2str(run)]);
else
    ppparams.resultfolder = ['SPMMAT-' task '_' params.SPMMAT_analysisname];
    ppparams.fmriresultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.SPMMAT_analysisname]);
end

if ~isfolder(ppparams.fmriresultmap)
    e.message = ['No SPMMAT forlder for subject ' num2str(sub) ' session ' num2str(ses)];
    error(e)
    return
end

load(fullfile(ppparams.fmriresultmap,'SPM.mat'));

peakX = max(SPM.xX.X,[],'all');

mask = spm_read_vols(spm_vol(fullfile(ppparams.fmriresultmap,'mask.nii')));
vmask = find(mask>0);

Bmean_dat = spm_read_vols(spm_vol(fullfile(ppparams.fmriresultmap,SPM.Vbeta(end).fname)));
vbmask = find(Bmean_dat(vmask)>0);

nCon = numel(SPM.xCon);

for icon=1:nCon
    VC = spm_vol(fullfile(ppparams.fmriresultmap,SPM.xCon(icon).Vcon.fname));
    Con_dat = spm_read_vols(VC);

    psc_dat = zeros(VC.dim);
    psc_dat(vmask(vbmask)) = (Con_dat(vmask(vbmask)) ./ Bmean_dat(vmask(vbmask))) * peakX * 100;

    PVC = VC;
    rmfield(PVC,'pinfo');
    PVC.fname = fullfile(ppparams.fmriresultmap,['PSC_' SPM.xCon(icon).Vcon.fname]);
    PVC.descrip = 'my_spmbatch - Percent SIgnal Change';
    PVC.dt = [spm_type('float32'),spm_platform('bigend')];
    PVC.n = [1 1];
    PVC = spm_write_vol(PVC,psc_dat);

    clear Con_dat VC psc_dat PVC
end

clear mask vmask Bmean_dat