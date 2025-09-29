function out = my_spmbatch_run_vbmpreprocessing(sub,ses,datpath,paramsfile)

load(paramsfile)

global spmpath
spmpath = params.spm_path;

if params.onVSC, [datpath,params] = before_run_VSC(datpath,sub,ses,params); end

fprintf('\nStart analysing the data\n');

try
    %% preprocess anatomical scans
    [delfiles,keepfiles] = my_spmbatch_preprocess_anat(sub,ses,datpath,params);

    % Clean up unnecessary files
    cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,'anat',params.save_folder);
catch e
    fprintf('\nPP_Error\n');
    fprintf('\nThe error was: \n%s\n',e.message)
end

if params.onVSC, [datpath,params] = after_run_VSC(datpath,sub,ses,params); end

fprintf('\nPP_Completed\n');

out = 1;