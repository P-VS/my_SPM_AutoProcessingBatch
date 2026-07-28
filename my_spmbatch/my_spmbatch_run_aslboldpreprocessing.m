function out = my_spmbatch_run_aslboldpreprocessing(sub,ses,run,task,datpath,paramsfile)

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
    
    %% preprocess ASL-BOLD scans
    if params.preprocess_asl
        [delfiles,keepfiles] = my_spmbatch_fasl(sub,ses,run,task,datpath,params);

        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,'perf',params.perf_save_folder); 
    end

    %% denoise functional scans
    if params.do_denoising
        [delfiles,keepfiles] = my_spmbatch_denoise(sub,ses,run,task,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,params.func_save_folder,params.perf_save_folder);
    end
catch e
    fprintf(['\nError processing ' num2str(sub,['%0' num2str(params.sub_digits) 'd']) ' ses-' num2str(ses,'%03d') ' run-' num2str(run,'%02d') ' task-' task '\n']);

    nlogfname = fullfile(datpath,['error_fmri_preprocessing_' num2str(sub,['%0' num2str(params.sub_digits) 'd']) '_ses-' num2str(ses,'%03d') '_run-' num2str(run,'%02d') '_task-' task '.txt']);

    fid = fopen(nlogfname, 'w');
    fprintf(fid,['Error processing ' num2str(sub,['%0' num2str(params.sub_digits) 'd']) '_ses-' num2str(ses,'%03d') '_run-' num2str(run,'%02d') '_task-' task '\n\n']);
    fprintf(fid,'\nThe error was: \n%s\n',e.message);
    fprintf(fid,'\n');
    if isfield(e,'stack')
        for istack=1:numel(e.stack)
            fprintf(fid,'\nError in file %s',e.stack(istack).file);
            fprintf(fid,'\name %s',e.stack(istack).name);
            fprintf(fid,' line %s\n',num2str(e.stack(istack).line));
        end
    end
    fclose(fid);

    fprintf('\nThe error was: \n%s\n',e.message)
    if isfield(e,'stack')
        for istack=1:numel(e.stack)
            fprintf('\nError in file %s',e.stack(istack).file);
            fprintf('\name %s',e.stack(istack).name);
            fprintf(' line %s\n',num2str(e.stack(istack).line));
        end
    end
end

if params.onVSC, [datpath,params] = after_run_VSC(datpath,sub,ses,params); end

fprintf('\nCompleted\n');

out = 1;