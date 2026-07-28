function AutoSPMpercentSignalChange

%Script to do change the Beta en contrast maps into Percent Signal Change maps in SPM12
%
%The individual first level analyses should have run
%
%* IMPORTANT: !! Look at your data before starting any (pre)processing. Losing time in trying to process bad data makes no sense !!

%Script written by dr. Peter Van Schuerbeek (Radiology UZ Brussel)

%% Give path to SPM25

clear 'all'

params.spm_path = '/Users/petervanschuerbeek/Library/Mobile Documents/com~apple~CloudDocs/Matlab/spm25';

%% Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/ASLBOLD_Manon/data'; 

sublist = [1:31]; %﻿list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
params.sub_digits = 2; %if 2 the subject folder is sub-01, if 3 the subject folder is sub-001, ...

nsessions = [1]; %nsessions>0
 
params.task = {'PREcog'}; %text string that is in between task_ and _bold in your fNRI nifiti filename

%% In case of multiple runs in the same session exist
params.func.mruns = false; %true if run number is in filename
params.func.runs = [2]; %the index of the runs (in filenames run-(index))

params.SPMMAT_analysisname = 'MEICA-ASL_SPLINE';
params.modality = 'fasl'; %'fmri' or 'fasl'
params.isaslbold = true;

params.onVSC = false; % !!!Only true if using the VSC with a VUB account!!!
params.use_parallel = true; 
params.maxprocesses = 5; %Best not too high to avoid memory problems

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

my_spmbatch_start_PSCanalysis(sublist,nsessions,datpath,params);

cd(curdir)

if ~params.onVSC
    if isfile(new_spm_read_vols_file), movefile(new_spm_read_vols_file,old_spm_read_vols_file); end
end

fprintf('\nDone with processing the data\n')

if params.onVSC, exit; else clear 'all'; end

end