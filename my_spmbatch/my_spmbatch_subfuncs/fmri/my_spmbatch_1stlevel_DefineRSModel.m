function [fmri_spec,ppparams] = my_spmbatch_1stlevel_DefineRSModel(sub,ses,run,task,datpath,params,ppparams)

%% fMRI model specification

jsondat = fileread(ppparams.frun(1).funcjsonfile);
jsondat = jsondecode(jsondat);

tr = jsondat.RepetitionTime;

if isfield(jsondat,'SliceTiming')
    SliceTiming = jsondat.SliceTiming;
    nsl= ceil(numel(SliceTiming)/numel(find(SliceTiming==SliceTiming(1))));
else
    nsl = 1;
end

fmri_spec.dir = {ppparams.resultmap};
fmri_spec.timing.units = 'secs';
fmri_spec.timing.RT = tr;
fmri_spec.timing.fmri_t = 1; 
fmri_spec.timing.fmri_t0 = 1; 

for ir=1:numel(params.iruns)
    for ne=1:numel(params.func.echoes)
        nsess = (ir-1)*numel(params.func.echoes)+ne;
    
        fmri_spec.sess(nsess).scans = ppparams.ppfmridat{ir}.sess{ne}.func(:,1);

        fmri_spec.sess(nsess).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    
        fmri_spec.sess(nsess).multi = {''};
        fmri_spec.sess(nsess).regress = struct('name', {}, 'val', {});
        fmri_spec.sess(nsess).multi_reg = {''};
    
        if ~contains(ppparams.frun(1).func(1).funcfile,'f'), fmri_spec.sess(nsess).hpf = params.hpf;
        else
            Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(1).func(1).funcfile));
            fmri_spec.sess(nsess).hpf = tr * (numel(Vfunc)-1);
        end
    end
end

fmri_spec.fact = struct('name', {}, 'levels', {});
fmri_spec.bases.hrf.derivs = [0 0]; 
fmri_spec.volt = 1;
fmri_spec.global = 'None';

Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(1).func(1).funcfile));
nvols = min([numel(Vfunc),50]);
fdata = spm_read_vols(Vfunc(1:nvols));
mask = my_spmbatch_mask(fdata);

Vmask = Vfunc(1);
rmfield(Vmask,'pinfo');
Vmask.fname = fullfile(ppparams.preprocfmridir,['mask_' ppparams.frun(1).func(1).funcfile]);
Vmask.descrip = 'my_spmbatch - mask';
Vmask.dt = [spm_type('float32'),spm_platform('bigend')];
Vmask.n = [1 1];
Vmask = spm_write_vol(Vmask,mask);

ppparams.mask_file = Vmask.fname;

clear fdata mask Vfunc

fmri_spec.mthresh = 0.8; 
fmri_spec.cvi = 'none';

fmri_spec.mask = {Vmask.fname};
