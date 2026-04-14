function AutoPPIanalyis

%Script to do the auto 1st level fMRI processing in SPM25
%
%Preparation:
%* Organise the data in BIDS format
%    - datpath
%        -sub-##
%            -ses-00# (if your experiment contains multiple session per subject)
%                -anat: containes the anatomical data (3D T1)
%                -func: containes the fmri data
%            
%* Make sure all data is preprocessed
%    - preproc_anat
%    - preproc_func
    
%* IMPORTANT: !! Look at your preprocessed data before starting any analysis. It makes no sense to lose time in trying to process bad data !!

%% Give path to SPM25

clear 'all'

params.spm_path = '/Users/accurad/Library/Mobile Documents/com~apple~CloudDocs/Matlab/spm25';

%% Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/data';  %'/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/data';

sublist = [1]; %﻿list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
params.sub_digits = 2; %if 2 the subject folder is sub-01, if 3 the subject folder is sub-001, ...

nsessions = [1]; %nsessions>0
 
params.task = {'PREcog'}; %text string that is in between task_ and _bold in your fNRI nifiti filename

%% In case of multiple runs in the same session exist
params.func.mruns = true; %true if run number is in filename
params.func.runs = [1]; %the index of the runs (in filenames run-(index))

params.SPMMAT_analysisname = 'MEICA-BOLD_OPTHRF';
params.modality = 'fmri'; %'fmri' of 'fasl'
params.isaslbold = true;

params.preprocfmridir = 'preproc_meica_bold'; %directory with the preprocessed fMRI data

params.onVSC = false; % !!!Only true if using the VSC with a VUB account!!!
params.use_parallel = true; 
params.run_background = true;
params.maxprocesses = 2; %Best not too high to avoid memory problems
params.keeplogs = false;

%% PPI/BSC data parameters

params.PPI_analysisname = 'MF_9ROIs';
params.VOIfolder = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/PPI_Regios_MF';

params.doPPI = true;
params.doBSC = false;

%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%--------------------------------------------------------------------------------
 
params.analysis_type = 'GLM';

global spmpath
spmpath = params.spm_path;

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

if ~params.onVSC
    restoredefaultpath

    if exist(params.spm_path,'dir'), addpath(genpath(params.spm_path)); end
    if exist(params.my_spmbatch_path,'dir'), addpath(genpath(params.my_spmbatch_path)); end
    
    old_spm_read_vols_file=fullfile(spm('Dir'),'spm_read_vols.m');
    new_spm_read_vols_file=fullfile(spm('Dir'),'old_spm_read_vols.m');
    
    if isfile(old_spm_read_vols_file), movefile(old_spm_read_vols_file,new_spm_read_vols_file); end
end
       
fprintf('Start with processing the data\n')

warnstate = warning;
warning off;

spm('defaults', 'FMRI');

curdir = pwd;

my_spmbatch_start_ppianalysis(sublist,nsessions,datpath,params);

cd(curdir)

if ~params.onVSC
    if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end
end

fprintf('\nDone with processing the data\n')

if params.onVSC, exit; else clear 'all'; end

end