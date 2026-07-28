function out = my_spmbatch_run_rsfmri1stlevel(sub,ses,run,task,datpath,paramsfile)

load(paramsfile)

global spmpath
spmpath = params.spm_path;

if params.onVSC, [datpath,params] = before_run_VSC(datpath,sub,ses,params); end

try
    %% make batch
    params = my_spmbatch_rsfmri1stlevel_processing(sub,ses,run,task,datpath,params);

catch e
    fprintf(['\nError processing ' num2str(sub,['%0' num2str(params.sub_digits) 'd']) ' ses-' num2str(ses,'%03d') ' run-' num2str(run,'%02d') ' task-' task '\n']);

    nlogfname = fullfile(datpath,['error_fmri_CONN1stlevel_' num2str(sub,['%0' num2str(params.sub_digits) 'd']) '_ses-' num2str(ses,'%03d') '_run-' num2str(run,'%02d') '_task-' task '.txt']);

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