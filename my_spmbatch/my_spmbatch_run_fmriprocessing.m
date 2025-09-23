function out = my_spmbatch_run_fmriprocessing(sub,ses,run,task,datpath,paramsfile)

load(paramsfile)

try
    %% do smoothing first (if not done before)
    if params.do_smoothing
        params = my_spmbatch_fmriprocessing_do_smoothing(sub,ses,run,task,datpath,params);
    end

    %% make batch
    if contains(params.analysis_type,'GLM') 
        params = my_spmbatch_fmrilevel1processing(sub,ses,run,task,datpath,params);
    elseif contains(params.analysis_type,'ICA') 
        params = my_spmbatch_ica1stlevel(sub,ses,run,task,datpath,params);
    end
catch e
    fprintf('\nPP_Error\n');
    fprintf('\nThe error was: \n%s\n',e.message)
end

fprintf('\nPP_Completed\n');

out = 1;