function [datpath,params] = after_run_VSC(datpath,sub,ses,params)

fprintf('\nStart copying results\n');

if isfolder(fullfile(params.new_subpath,'preproc_anat'))
    copyfile(fullfile(params.new_subpath,'preproc_anat'),fullfile(params.orig_subpath,'preproc_anat')); 
    rmdir(fullfile(params.new_subpath,'preproc_anat'),'s');
end
if isfield(params,'func_save_folder')
    if isfolder(fullfile(params.new_subpath,params.func_save_folder))
        copyfile(fullfile(params.new_subpath,params.func_save_folder),fullfile(params.orig_subpath,params.func_save_folder)); 
        rmdir(fullfile(params.new_subpath,params.func_save_folder),'s');
    end
end
if isfield(params,'perf_save_folder')
    if isfolder(fullfile(params.new_subpath,params.perf_save_folder))
        copyfile(fullfile(params.new_subpath,params.perf_save_folder),fullfile(params.orig_subpath,params.perf_save_folder)); 
        rmdir(fullfile(params.new_subpath,params.perf_save_folder),'s');
    end
end
if isfolder(fullfile(params.new_subpath,'anat')), rmdir(fullfile(params.new_subpath,'anat'),'s'); end
if isfolder(fullfile(params.new_subpath,'fmap')), rmdir(fullfile(params.new_subpath,'fmap'),'s'); end
if isfolder(fullfile(params.new_subpath,'func')), rmdir(fullfile(params.new_subpath,'func'),'s'); end
if isfolder(fullfile(params.new_subpath,'perf')), rmdir(fullfile(params.new_subpath,'perf'),'s'); end

datpath = params.orig_subpath;