function [ppparams,delfiles,keepfiles] = my_spmbatch_aslbold_normalization(ppparams,params,delfiles,keepfiles)

Vcbf = spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).cbfprefix ppparams.perf(1).cbffile]));

for i=1:numel(Vcbf)
    wcbffiles{i,1} = [Vcbf(i).fname ',' num2str(i)];
end

%% Normalization of the M0 scan

ppparams.deffile = fullfile(ppparams.subperfdir,['y_' ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile ',1']);

m0normest.subj.vol = {fullfile(ppparams.subperfdir,[ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile ',1'])};
m0normest.eoptions.biasreg = 0.0001;
m0normest.eoptions.biasfwhm = 60;
m0normest.eoptions.tpm = {fullfile(params.spm_path,'tpm','TPM.nii')};
m0normest.eoptions.affreg = 'mni';
m0normest.eoptions.reg = [0 0 0.1 0.01 0.04];
m0normest.eoptions.fwhm = 0;
m0normest.eoptions.samp = 3;

spm_run_norm(m0normest);

delfiles{numel(delfiles)+1} = {ppparams.deffile};

%Write the spatially normalised  m0scan data

m0normw.woptions = spm_get_defaults('normalise.write');
m0normw.woptions.vox = params.func.normvox;

dt = Vcbf(1).dt;
if dt(1)==spm_type('uint16')
    m0normw.woptions.interp = 4;
else
    m0normw.woptions.interp = 1;
end

m0normw.subj.def = {ppparams.deffile};
m0normw.subj.resample = {fullfile(ppparams.subperfdir,[ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile ',1'])};

spm_run_norm(m0normw);

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subperfdir,['w' ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile])};
ppparams.perf(1).wm0file = ['w' ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile];


%% Normalise CBF data

%Write the spatially normalised  CBF data

cbfnormw = m0normw;
cbfnormw.subj.def = {ppparams.deffile};
cbfnormw.subj.resample = wcbffiles(:,1);

spm_run_norm(cbfnormw);

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subperfdir,['w' ppparams.perf(1).cbfprefix ppparams.perf(1).cbffile])};

ppparams.perf(1).wcbffile = ['w' ppparams.perf(1).cbfprefix ppparams.perf(1).cbffile];

clear Vcbf

%Write the spatially normalised  mean CBF data

Vmcbf = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).meancbf));

wmcbffiles{1,1} = [Vmcbf(1).fname ',1'];

mcbfnormw = cbfnormw;
mcbfnormw.subj.resample = wmcbffiles(:,1);

spm_run_norm(mcbfnormw);

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subperfdir,['w' ppparams.perf(1).meancbf])};

ppparams.perf(1).wmcbffile = ['w' ppparams.perf(1).meancbf];

clear Vmcbf