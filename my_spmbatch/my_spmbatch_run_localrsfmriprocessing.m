function out = my_spmbatch_run_localrsfmriprocessing(sub,ses,run,task,datpath,paramsfile)

load(paramsfile)

global spmpath
spmpath = params.spm_path;

if params.onVSC, [datpath,params] = before_run_VSC(datpath,sub,ses,params); end

try
    %% make batch
    params = my_spmbatch_localrsfmriprocessing(sub,ses,run,task,datpath,params);
catch e
    fprintf('\nPP_Error\n');
    fprintf('\nThe error was: \n%s\n',e.message)
end

if params.onVSC, [datpath,params] = after_run_VSC(datpath,sub,ses,params); end

fprintf('\nPP_Completed\n');

out = 1;