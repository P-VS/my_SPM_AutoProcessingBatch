function [ppparams,delfiles,keepfiles] = my_spmbatch_dosmoothasl(ppparams,params,delfiles,keepfiles)

%% Smooth CBF series
Vcbf = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).wcbffile));
cbfdat = spm_read_vols(Vcbf);

tdim = numel(Vcbf);

Vout = Vcbf;
for iv=1:tdim
    scbfdat = my_spmbatch_smooth(cbfdat(:,:,:,iv),Vcbf(iv),[],[params.func.smoothfwhm params.func.smoothfwhm params.func.smoothfwhm],0);

    Vout(iv).fname = fullfile(ppparams.subperfdir,['s' ppparams.perf(1).wcbffile]);
    Vout(iv).descrip = 'my_spmbatch - smooth';
    Vout(iv).n = [iv 1];
    Vout(iv) = spm_write_vol(Vout(iv),scbfdat);

    if mod(iv,50)==0, fprintf(['done smoothing vols ' num2str(iv) ' of ' num2str(tdim) '\n']); end
end

clear Vout cbfdat

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subperfdir,['s' ppparams.perf(1).wcbffile])};    

ppparams.perf(1).scbffile = ['s' ppparams.perf(1).wcbffile];

clear Vcbf

%% Smooth mean CBF

Vmcbf = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).wmcbffile));
mcbfdat = spm_read_vols(Vmcbf);

smcbfdat = my_spmbatch_smooth(mcbfdat,Vmcbf,[],[params.func.smoothfwhm params.func.smoothfwhm params.func.smoothfwhm],0);

Vout = Vmcbf;
Vout.fname = fullfile(ppparams.subperfdir,['s' ppparams.perf(1).wmcbffile]);
Vout.descrip = 'my_spmbatch - smooth';
Vout.n = [1 1];
Vout = spm_write_vol(Vout,smcbfdat);

clear sasldat Vout asldat Vmcbf

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subperfdir,['s' ppparams.perf(1).wmcbffile])};    

ppparams.perf(1).swmcbffile = ['s' ppparams.perf(1).wmcbffile];
