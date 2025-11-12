function AutoLocalRSfMRI_fMRI

%Script to do the auto local RS fMRI processing in REST ((f)ALFF and ReHo)
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

datpath = '/Volumes/LaCie/UZ_Brussel/ME_fMRI_GE/data';  %'/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/data';

sublist = [1,2,4:11]; %﻿list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
params.sub_digits = 2; %if 2 the subject folder is sub-01, if 3 the subject folder is sub-001, ...

nsessions = [1]; %nsessions>0
 
params.task = {'ME-EFT'}; %text string that is in between task_ and _bold in your fNRI nifiti filename

params.analysisname = 'ME-DUNE-NonBOLD';
params.modality = 'fmri'; %'fmri' of 'fasl'
params.isaslbold = false;

params.onVSC = false; % !!!Only true if using the VSC with a VUB account!!!
params.use_parallel = false; 
params.run_background = true;
params.maxprocesses = 2; %Best not too high to avoid memory problems
params.keeplogs = false;

%% fMRI data parameters
    params.preprocfmridir = 'preproc_func_ME-DUNE-NonBOLD'; %directory with the preprocessed fMRI data
    params.fmri_prefix = 'wcdavure'; %fMRI file name of form [fmri_prefix 'sub-ii_task-..._' fmri_endfix '.nii']
    
    %In case of multiple runs in the same session exist
    params.func.mruns = false; %true if run number is in filename
    params.func.runs = [1]; %the index of the runs (in filenames run-(index))
    params.func.use_runs = 'separately'; % 'separately' or 'together' (this parameter is ignored if mruns is false)
    %'separately': a separate analysis is done per run (default=separately)
    %'together': all runs are combined in 1 analysis 
    
    % For ME-fMRI
    params.func.meepi = true; %true if echo number is in filename
    params.func.echoes = [1]; %the index of echoes in ME-fMRI used in the analysis. If meepi=false, echoes=[1]. 

%% Analysis parameters
    params.do_ALFF = true;
    params.do_fALFF = true;
    params.do_ReHo = true; %KCC Kendall’s coefficient of concordance for N neighbouring voxels

    params.LF_band = [0.01 0.1]; %BOLD frequency band for (f)ALFF (default=[0.008 0.1])
    params.Nvoxels = 27; %number of neighbouring voxels (27, 19, or 7) (default=27)

%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%--------------------------------------------------------------------------------
 
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

my_spmbatch_start_LocalRSfmriprocessing(sublist,nsessions,datpath,params);

cd(curdir)

if ~params.onVSC
    if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end
end

fprintf('\nDone with processing the data\n')

if params.onVSC, exit; else clear 'all'; end

end