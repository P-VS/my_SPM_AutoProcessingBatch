function make_AALROImasks

[my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));
fsplit = split(my_spmbatch_path,'small_scripts');

AALim_file = fullfile(fsplit{1},'templates','AAL3v1_1mm.nii');
AALim_txt = fullfile(fsplit{1},'templates','AAL3v1_1mm.txt');

VAAL = spm_vol(AALim_file);
AALim_dat = spm_read_vols(VAAL);

AALlist = importdata(AALim_txt);

for iroi=1:size(AALlist.textdata,1)
    tmp = find(and(AALim_dat>iroi-0.5,AALim_dat<iroi+0.5));
    roimask = zeros(VAAL.dim(1),VAAL.dim(2),VAAL.dim(3));
    roimask(tmp) = 1;

    AALR = VAAL;
    AALR.fname =  fullfile(fsplit{1},'templates','AALRois',[num2str(str2num(AALlist.textdata{iroi,1}),'%03d') '_' AALlist.textdata{iroi,2} '.nii']);
    AALR=spm_write_vol(AALR,roimask);

    clear AALR roimask tmp
end

mask = zeros(VAAL.dim(1),VAAL.dim(2),VAAL.dim(3));
mask(AALim_dat>0.1) = 1;

VAAL.fname =  fullfile(fsplit{1},'templates','AALRois','AALmask.nii');
VAAL=spm_write_vol(VAAL,mask);

clear VAAL mask
