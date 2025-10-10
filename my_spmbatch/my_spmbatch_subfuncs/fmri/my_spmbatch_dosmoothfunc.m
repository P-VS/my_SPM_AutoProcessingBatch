function [ppparams,delfiles,keepfiles] = my_spmbatch_dosmoothfunc(ppparams,params,ie,delfiles,keepfiles)
 
Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
funcdat = spm_read_vols(Vfunc);

tdim = numel(Vfunc);

Vout = Vfunc;
for iv=1:tdim
    sfuncdat = my_spmbatch_smooth(funcdat(:,:,:,iv),Vfunc(iv),[],[params.func.smoothfwhm params.func.smoothfwhm params.func.smoothfwhm],0);

    Vout(iv).fname = fullfile(ppparams.subfuncdir,['s' ppparams.func(ie).prefix ppparams.func(ie).funcfile]);
    Vout(iv).descrip = 'my_spmbatch - smooth';
    Vout(iv).n = [iv 1];
    Vout(iv) = spm_write_vol(Vout(iv),sfuncdat);

    clear sfuncdat
end

clear Vout funcdat

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subfuncdir,['s' ppparams.func(ie).prefix ppparams.func(ie).funcfile])};    

ppparams.func(ie).prefix = ['s' ppparams.func(ie).prefix];