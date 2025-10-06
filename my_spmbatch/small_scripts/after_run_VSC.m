function [datpath,params] = after_run_VSC(datpath,sub,ses,params)

fprintf('\nStart copying results\n');

if isfolder(fullfile(params.new_subpath,'preproc_anat')), copyfile(fullfile(params.new_subpath,'preproc_anat'),fullfile(params.orig_subpath,'preproc_anat')); end
if isfield(params,'func_save_folder')
    if isfolder(fullfile(params.new_subpath,params.func_save_folder)), copyfile(fullfile(params.new_subpath,params.func_save_folder),fullfile(params.orig_subpath,params.func_save_folder)); end
end
if isfield(params,'perf_save_folder')
    if isfolder(fullfile(params.new_subpath,params.perf_save_folder)), copyfile(fullfile(params.new_subpath,params.perf_save_folder),fullfile(params.orig_subpath,params.perf_save_folder)); end
end

datpath = params.orig_subpath;