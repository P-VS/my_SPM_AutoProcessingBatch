function out = my_spmbatch_run_aslpreprocessing(sub,ses,run,datpath,paramsfile)

load(paramsfile)

try
    %% preprocess anatomical scans
    if params.preprocess_anatomical
        [delfiles,keepfiles] = my_spmbatch_preprocess_anat(sub,ses,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,'anat','preproc_anat');
    end
    
    %% preprocess asl scans
    if params.preprocess_pcasl
        [delfiles,keepfiles] = my_spmbatch_aslpreprocessed(sub,ses,run,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,'perf',params.save_folder);
    end
catch e
    fprintf(['\nError processing ' num2str(sub,['%0' num2str(params.sub_digits) 'd']) ' ses-' num2str(ses,'%03d') ' run-' num2str(run,'%02d') ' task-' task '\n']);

    nlogfname = fullfile(datpath,['error_asl_preprocessing_' num2str(sub,['%0' num2str(params.sub_digits) 'd']) '_ses-' num2str(ses,'%03d') '_run-' num2str(run,'%02d') '_task-' task '.txt']);

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

fprintf('\nCompleted\n');

out = 1;