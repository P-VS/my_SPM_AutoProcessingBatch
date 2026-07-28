function AutoCONN_1stlevel_RSFMRI

%Script to do the auto 1st level ROI-ROI and seed-voxel analysis for RS-fMRI in CON
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

params.spm_path = '/Users/petervanschuerbeek/Library/Mobile Documents/com~apple~CloudDocs/Matlab/spm25';
params.conn_path = '/Users/petervanschuerbeek/Library/Mobile Documents/com~apple~CloudDocs/Matlab/conn';

%% Give the basic input information of your data

datpath = '/Users/petervanschuerbeek/fMRI_data/HumanIT/IndividueleData'; 

sublist = [2,3]; %﻿list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
params.sub_digits = 2; %if 2 the subject folder is sub-01, if 3 the subject folder is sub-001, ...

nsessions = [1,2]; %nsessions>0

params.task = {'rest'}; %text string that is in between task_ and _bold in your fNRI nifiti filename

params.analysisname = 'RS-fMRI';

params.onVSC = false; % !!!Only true if using the VSC with a VUB account!!!
params.use_parallel = false; 
params.maxprocesses = 3; %Best not too high to avoid memory problems

%% fMRI data parameters
params.preprocfmridir = 'preproc_func'; %directory with the preprocessed fMRI data
params.fmri_prefix = 'swdfavure'; %fMRI file name of form [fmri_prefix 'sub-ii_task-..._' fmri_endfix '.nii']

%In case of multiple runs in the same session exist
params.func.mruns = false; %true if run number is in filename
params.func.runs = [1]; %the index of the runs (in filenad mes run-(index))

% For ME-fMRI
params.func.meepi = false; %true if echo number is in filename
params.func.echoes = [1:3]; %the index of echoes in ME-fMRI used in the analysis. If meepi=false, echoes=[1]. 

%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%--------------------------------------------------------------------------------

global spmpath
spmpath = params.spm_path;

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

if ~params.onVSC
    restoredefaultpath

    if exist(params.spm_path,'dir'), addpath(genpath(params.spm_path)); end
    if exist(params.conn_path,'dir'), addpath(genpath(params.conn_path)); end
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

my_spmbatch_start_1stlevel_RSFMRI(sublist,nsessions,datpath,params);

cd(curdir)

if ~params.onVSC
    if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end
end

fprintf('\nDone with processing the data\n')

if params.onVSC, exit; else clear 'all'; end

end