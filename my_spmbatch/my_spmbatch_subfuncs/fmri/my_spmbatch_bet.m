function [outfile,delfiles] = my_spmbatch_bet(infolder,infile,ppparams,params,delfiles,keepfiles)

% ---------------------------------------------------------------------
% Brain extraction using MRTOOL
% ---------------------------------------------------------------------
fprintf('Start brain extraction \n')

fprintf('Do segmentation \n')

preproc.channel.vols = {fullfile(infolder,[infile ',1'])};
preproc.channel.biasreg = 0.001;
preproc.channel.biasfwhm = 60;
preproc.channel.write = [0 0];
preproc.tissue(1).tpm = {fullfile(params.spm_path,'tpm','TPM.nii,1')};
preproc.tissue(1).ngaus = 1;
preproc.tissue(1).native = [1 0];
preproc.tissue(1).warped = [0 0];
preproc.tissue(2).tpm = {fullfile(params.spm_path,'tpm','TPM.nii,2')};
preproc.tissue(2).ngaus = 1;
preproc.tissue(2).native = [1 0];
preproc.tissue(2).warped = [0 0];
preproc.tissue(3).tpm = {fullfile(params.spm_path,'tpm','TPM.nii,3')};
preproc.tissue(3).ngaus = 2;
preproc.tissue(3).native = [1 0];
preproc.tissue(3).warped = [0 0];
preproc.tissue(4).tpm = {fullfile(params.spm_path,'tpm','TPM.nii,4')};
preproc.tissue(4).ngaus = 3;
preproc.tissue(4).native = [0 0];
preproc.tissue(4).warped = [0 0];
preproc.tissue(5).tpm = {fullfile(params.spm_path,'tpm','TPM.nii,5')};
preproc.tissue(5).ngaus = 4;
preproc.tissue(5).native = [0 0];
preproc.tissue(5).warped = [0 0];
preproc.tissue(6).tpm = {fullfile(params.spm_path,'tpm','TPM.nii,6')};
preproc.tissue(6).ngaus = 2;
preproc.tissue(6).native = [0 0];
preproc.tissue(6).warped = [0 0];
preproc.warp.mrf = 1;
preproc.warp.cleanup = 1;
preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
preproc.warp.affreg = 'mni';
preproc.warp.fwhm = 0;
preproc.warp.samp = 3;
preproc.warp.write = [0 1];
preproc.warp.vox = NaN;
preproc.warp.bb = [NaN NaN NaN;NaN NaN NaN];

spm_preproc_run(preproc);

fc1im = fullfile(infolder,['c1' infile]);
fc2im = fullfile(infolder,['c2' infile]);
fc3im = fullfile(infolder,['c3' infile]);

sname = split(infile,'.nii');
delfiles{numel(delfiles)+1} = {fullfile(infolder,[sname{1} '._seg8.mat'])};
delfiles{numel(delfiles)+1} = {fullfile(infolder,['y_' sname{1} '.nii'])};  

gm = spm_read_vols(spm_vol(fc1im));
wm = spm_read_vols(spm_vol(fc2im));
csf = spm_read_vols(spm_vol(fc3im));

Vin = spm_vol(fullfile(infolder,infile));
indat = spm_read_vols(Vin(1));

indat((gm+wm+csf)<0.1) = 0;

outfile = ['b' infile];

Vout = Vin(1);
Vout.fname = fullfile(infolder,outfile);
Vout = spm_write_vol(Vout,indat);

clear Vout indat gm wm csf

delfiles{numel(delfiles)+1} = {fullfile(infolder,outfile)};

fprintf('Done brain extraction \n')