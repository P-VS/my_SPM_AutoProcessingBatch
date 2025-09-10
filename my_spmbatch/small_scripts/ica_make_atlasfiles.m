function ica_make_atlasfiles

output_directory = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/AAtlasROIS';

network(1).roinumbers = [80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90];
network(1).label = 'Visuospatial_Network';

network(2).roinumbers = [64, 65, 66, 67, 68, 69];
network(2).label = 'Sensorimotor_Network';

network(3).roinumbers = [1, 2, 3, 4, 5, 6, 7, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51];
network(3).label = 'Salience_Network';

network(4).roinumbers = [34, 35, 36, 37, 38, 39, 58, 59, 60, 61, 62, 63];
network(4).label = 'ExecutiveControl_Network';

network(5).roinumbers = [27, 28, 29, 30, 31, 32, 33];
network(5).label = 'Language_Network';

network(6).roinumbers = [16, 17, 18, 19, 20, 21, 22, 23, 24, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79];
network(6).label = 'DefaultMode_Network';

network(7).roinumbers = [1:90];
network(7).label = 'All';

[my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));
fsplit = split(my_spmbatch_path,'small_scripts');
template_file = fullfile(fsplit{1},'templates','RSN.nii');

Vtemp = spm_vol(template_file);
templateIm = spm_read_vols(Vtemp);

for inet=1:numel(network)
    netwim = sum(templateIm(:,:,:,network(inet).roinumbers),4);
    netwim(netwim>0) = 1;

    Vout = Vtemp(1);
    Vout.fname=fullfile(output_directory,['network-' num2str(inet,'%03d') '_' network(inet).label '.nii']);
    Vout=spm_write_vol(Vout,netwim);

    clear Vout netwim
end