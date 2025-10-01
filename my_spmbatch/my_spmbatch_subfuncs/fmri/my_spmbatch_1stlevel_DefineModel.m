function [fmri_spec,ppparams] = my_spmbatch_1stlevel_DefineModel(sub,ses,run,task,datpath,params,ppparams)

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
fmri_spec.timing.fmri_t = nsl;
fmri_spec.timing.fmri_t0 = 1;

for ir=1:numel(params.iruns)
    % correct events file for dummy scans if needed
    dummys = floor(params.dummytime/tr);
    numparams = 0;

    T = readtable(ppparams.frun(ir).functsvfile,'FileType','text');
    ppparams.edat{ir}.onset = T.onset;
    ppparams.edat{ir}.duration = T.duration;
    ppparams.edat{ir}.trial_type = T.trial_type;
    if params.add_parametricModulation
        fnames = fieldnames(T);
        if numel(fnames)>3
            for ifield=4:numel(fnames)-3
                ppparams.edat{ir}.weight{ifield-3}.name = fnames{ifield};
                ppparams.edat{ir}.weight{ifield-3}.value = getfield(T,fnames{ifield});
            end
            numparams = numel(fnames)-6;
        else params.add_parametricModulation=false; end
    else params.add_parametricModulation=false;
    end

    ppparams.edat{ir}.onset = ppparams.edat{ir}.onset-dummys*tr;
    if params.isaslbold, ppparams.edat{ir}.onset=ppparams.edat{ir}.onset-(params.asl.LabelingDuration+params.asl.PostLabelDelay); end
    
    for it=1:numel(ppparams.edat{ir}.trial_type(:,1))
        ntrial_type(it,1) = convertCharsToStrings(ppparams.edat{ir}.trial_type(it,:));
    end
    
    ppparams.edat{ir}.trial_type = ntrial_type;
    
    [~,edatorder] = sort(ppparams.edat{ir}.trial_type);
    ppparams.edat{ir}.onset = ppparams.edat{ir}.onset(edatorder);
    ppparams.edat{ir}.duration = ppparams.edat{ir}.duration(edatorder);
    ppparams.edat{ir}.trial_type = ppparams.edat{ir}.trial_type(edatorder,:);
    if isfield(ppparams.edat{ir},'weight')
        for iw=1:numel(ppparams.edat{ir}.weight)
            ppparams.edat{ir}.weight{iw}.value = ppparams.edat{ir}.weight{iw}.value(edatorder); 
        end
    end
    
    numc=0;
    for trial=1:numel(ppparams.edat{ir}.onset)
        if isstring(ppparams.edat{ir}.trial_type(trial,:))
            trial_type = ppparams.edat{ir}.trial_type(trial,:);
        else
            trial_type = num2str(ppparams.edat{ir}.trial_type(trial,:));
        end
    
        if numc>0
            numc = numel(ppparams.edat{ir}.conditions)+1;
            for nc=1:numel(ppparams.edat{ir}.conditions)
                if strcmp(ppparams.edat{ir}.conditions{nc}.name,strtrim(trial_type)); numc=nc; end
            end
            if numc<numel(ppparams.edat{ir}.conditions)+1
                ppparams.edat{ir}.conditions{numc}.onsets = [ppparams.edat{ir}.conditions{numc}.onsets ppparams.edat{ir}.onset(trial)];
                ppparams.edat{ir}.conditions{numc}.durations = [ppparams.edat{ir}.conditions{numc}.durations ppparams.edat{ir}.duration(trial)];
                if params.add_parametricModulation
                    ppparams.edat{ir}.conditions{numc}.weight{ifield}.name = ppparams.edat{ir}.weight{ifield}.name; 
                    ppparams.edat{ir}.conditions{numc}.weight{ifield}.values = [ppparams.edat{ir}.conditions{numc}.weight{ifield}.values ppparams.edat{ir}.weight{ifield}.value(trial)]; 
                end
            else
                ppparams.edat{ir}.conditions{numc}.name = strtrim(trial_type);
                ppparams.edat{ir}.conditions{numc}.onsets = [ppparams.edat{ir}.onset(trial)];
                ppparams.edat{ir}.conditions{numc}.durations = [ppparams.edat{ir}.duration(trial)];
                if params.add_parametricModulation
                    for ifield=1:numparams
                        ppparams.edat{ir}.conditions{numc}.weight{ifield}.name = ppparams.edat{ir}.weight{ifield}.name; 
                        ppparams.edat{ir}.conditions{numc}.weight{ifield}.values = [ppparams.edat{ir}.weight{ifield}.value(trial)]; 
                    end
                end
            end
        else
            ppparams.edat{ir}.conditions{1}.name = strtrim(trial_type);
            ppparams.edat{ir}.conditions{1}.onsets = ppparams.edat{ir}.onset(trial);
            ppparams.edat{ir}.conditions{1}.durations = ppparams.edat{ir}.duration(trial);
            if params.add_parametricModulation
                for ifield=1:numparams
                    ppparams.edat{ir}.conditions{1}.weight{ifield}.name = ppparams.edat{ir}.weight{ifield}.name; 
                    ppparams.edat{ir}.conditions{1}.weight{ifield}.values = ppparams.edat{ir}.weight{ifield}.value(trial); 
                end
            end
            numc=1;
        end
    end
    
    for ne=1:numel(params.func.echoes)
        nsess = (ir-1)*numel(params.func.echoes)+ne;
    
        fmri_spec.sess(nsess).scans = ppparams.ppfmridat{ir}.sess{ne}.func(:,1);
    
        for nc=1:numel(ppparams.edat{1}.conditions)
            fmri_spec.sess(nsess).cond(nc).name = char(ppparams.edat{ir}.conditions{nc}.name);
            fmri_spec.sess(nsess).cond(nc).onset = ppparams.edat{ir}.conditions{nc}.onsets;
            fmri_spec.sess(nsess).cond(nc).duration = ppparams.edat{ir}.conditions{nc}.durations;
            fmri_spec.sess(nsess).cond(nc).tmod = 0;
            if ~params.add_parametricModulation
                fmri_spec.sess(nsess).cond(nc).pmod = struct('name', {}, 'param', {}, 'poly', {});
            else
                for ifield=1:numparams
                    fmri_spec.sess(nsess).cond(nc).pmod(ifield).name = ppparams.edat{ir}.conditions{nc}.weight{ifield}.name;
                    fmri_spec.sess(nsess).cond(nc).pmod(ifield).param = ppparams.edat{ir}.conditions{nc}.weight{ifield}.values;
                    fmri_spec.sess(nsess).cond(nc).pmod(ifield).poly = 1;
                end
            end
            fmri_spec.sess(nsess).cond(nc).orth = 1;
        end
    
        fmri_spec.sess(nsess).multi = {''};

        if params.isaslbold %contains(params.modality,'fasl') && contains(params.whichfile,'asl')
            labels = zeros(1,numel(ppparams.ppfmridat{ir}.sess{ne}.func));
            labels(2:2:end) = 1;
            fmri_spec.sess(nsess).regress.name = 'labeling';
            fmri_spec.sess(nsess).regress.val = labels;
        else
            fmri_spec.sess(nsess).regress = struct('name', {}, 'val', {});
        end
    
        fmri_spec.sess(nsess).multi_reg = {ppparams.frun(ir).confoundsfile};

        if ~contains(ppparams.frun(1).func(1).funcfile,'f'), fmri_spec.sess(nsess).hpf = params.hpf;
        else
            Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(1).func(1).funcfile));
            fmri_spec.sess(nsess).hpf = tr * (numel(Vfunc)-1);
        end
    end
end

fmri_spec.fact = struct('name', {}, 'levels', {});
if ~params.add_derivatives; fmri_spec.bases.hrf.derivs = [0 0]; 
else fmri_spec.bases.hrf.derivs = [1 1]; end
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

clear fdata mask

if contains(params.modality,'fasl') && contains(params.whichfile,'cbf')
    fmri_spec.mthresh = -1.0;
else
    fmri_spec.mthresh = 0.8; 
end
if contains(params.modality,'fasl')
    fmri_spec.cvi = params.model_serial_correlations;
else
    fmri_spec.cvi = 'AR(1)';
end
fmri_spec.mask = {Vmask.fname};
