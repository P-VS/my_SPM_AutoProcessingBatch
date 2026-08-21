function AutoCONN_groupSBC_RSFMRI

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

datpath = '/Volumes/LaCie/UZ_Brussel/HumanIT Kevin-Elke/data'; 

sublist = [2,3]; %﻿list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
params.sub_digits = 2; %if 2 the subject folder is sub-01, if 3 the subject folder is sub-001, ...

nsessions = [1,2]; %nsessions>0

params.task = {'rest'}; %text string that is in between task_ and _bold in your fNRI nifti filename

params.outfolder = '/Volumes/LaCie/UZ_Brussel/HumanIT Kevin-Elke/Group_CONN'; 
params.analysisname = 'CONN-HumanIT-2Subjects';

params.use_parallel = false; 
params.maxprocesses = 2; %Best not too high to avoid memory problems

%% fMRI data parameters
params.preprocfmridir = 'preproc_func'; %directory with the preprocessed fMRI data
params.fmri_prefix = 'swdfavure'; %fMRI file name of form [fmri_prefix 'sub-ii_task-..._' fmri_endfix '.nii']

%% In case of multiple runs in the same session exist
params.func.mruns = false; %true if run number is in filename
params.func.runs = [1]; %the index of the runs (in filenames run-(index))

%% CONN parameters
params.add = false; %add new data to an existing analysis

params.do_ROI2ROI = true;
params.do_Seed2voxel = true;
params.do_voxel2voxel = true;

%defines 2nd-level covariates (arbitrary continuous/categorical/ordinal data for each subject) 
params.sub_covariates{1}.names = 'All subjects';
params.sub_covariates{1}.vector = ones(numel(sublist),1); %subjects.effects{neffect} vector of size [nsubjects,1] defining second-level effects

%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%--------------------------------------------------------------------------------

global spmpath
spmpath = params.spm_path;

[params.my_spmbatch_path,~,~] = fileparts(mfilename('fullpath'));

restoredefaultpath

if exist(params.spm_path,'dir'), addpath(genpath(params.spm_path)); end
if exist(params.conn_path,'dir'), addpath(genpath(params.conn_path)); end
if exist(params.my_spmbatch_path,'dir'), addpath(genpath(params.my_spmbatch_path)); end

old_spm_read_vols_file=fullfile(spm('Dir'),'spm_read_vols.m');
new_spm_read_vols_file=fullfile(spm('Dir'),'old_spm_read_vols.m');

if isfile(old_spm_read_vols_file), movefile(old_spm_read_vols_file,new_spm_read_vols_file); end

fprintf('Start with processing the data\n')

warnstate = warning;
warning off;

spm('defaults', 'FMRI');

curdir = pwd;

my_spmbatch_conn_GroupSBC(sublist,nsessions,datpath,params);

cd(curdir)

if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end

fprintf('\nDone with processing the data\n')

clear 'all'

end