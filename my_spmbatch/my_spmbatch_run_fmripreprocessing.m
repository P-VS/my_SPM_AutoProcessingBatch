function out = my_spmbatch_run_fmripreprocessing(sub,ses,run,task,datpath,paramsfile)

load(paramsfile)

global spmpath
spmpath = params.spm_path;

if params.onVSC, [datpath,params] = before_run_VSC(datpath,sub,ses,params); end

fprintf('\nStart analysing the data\n');

try
    %% preprocess anatomical scans
    if params.preprocess_anatomical
        [delfiles,keepfiles] = my_spmbatch_preprocess_anat(sub,ses,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,'anat','preproc_anat');
    end
    
    %% preprocess functional scans
    if params.preprocess_functional
        [delfiles,keepfiles] = my_spmbatch_functional(sub,ses,run,task,datpath,params);
        
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,'func',params.func_save_folder);  
    end
    
    %% denoise functional scans
    if params.do_denoising && ~(params.preprocess_functional && params.func.denoise)
        [delfiles,keepfiles] = my_spmbatch_denoise(sub,ses,run,task,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,params.func_save_folder,params.func_save_folder);
    end
catch e
    fprintf('\nPP_Error\n');
    fprintf('\nThe error was: \n%s\n',e.message)
end

if params.onVSC, [datpath,params] = after_run_VSC(datpath,sub,ses,params); end

fprintf('\nPP_Completed\n');

out = 1;