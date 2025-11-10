function [ppparams,delfiles,keepfiles] = my_spmbatch_oldnormalization(ne,ppparams,params,delfiles,keepfiles)

Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ne).prefix ppparams.func(ne).funcfile]));

for i=1:numel(Vfunc)
    wfuncfiles{i,1} = [Vfunc(i).fname ',' num2str(i)];
end

%% Normalization of the functional scan
reffile = fullfile(ppparams.subfuncdir,ppparams.reffunc);
ppparams.deffile = spm_file(reffile, 'suffix','_sn', 'ext','.mat');
if ne==ppparams.echoes(1)
    if ~exist(ppparams.deffile,"file")
        template = {fullfile(params.spm_path,'toolbox','OldNorm','EPI.nii')};

        funcnormest.subj.source = {reffile};
        funcnormest.subj.wtsrc = '';
        funcnormest.eoptions.weight = '';
        funcnormest.eoptions.smosrc = 8;
        funcnormest.eoptions.smoref = 0;
        funcnormest.eoptions.regtype = 'mni';
        funcnormest.eoptions.reg = 1;
        funcnormest.eoptions.cutoff = 25;
        funcnormest.eoptions.nits = 16;
    
        spm_normalise(char(template),char(funcnormest.subj.source), ppparams.deffile,...
            char(funcnormest.eoptions.weight), char(funcnormest.subj.wtsrc), funcnormest.eoptions);
    end

    delfiles{numel(delfiles)+1} = {ppparams.deffile};
end

%% Normalise func

%Write the spatially normalised data

funcnormw.woptions.vox = params.func.normvox;
funcnormw.woptions.preserve = 0;
funcnormw.woptions.bb = [[-78 -112  -70]; [78   76   86]];
funcnormw.woptions.wrap = [0 0 0];
funcnormw.woptions.prefix = 'w';

dt = Vfunc(1).dt;
if ~(dt(1)==spm_type('uint16'))
    funcnormw.woptions.interp = 1;
end

funcnormw.subj.resample = wfuncfiles(:,1);

spm_write_sn(char(funcnormw.subj.resample), ppparams.deffile, funcnormw.woptions);

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subfuncdir,['w' ppparams.func(ne).prefix ppparams.func(ne).funcfile])};

ppparams.func(ne).prefix = ['w' ppparams.func(ne).prefix];

clear Vfunc Vtemp