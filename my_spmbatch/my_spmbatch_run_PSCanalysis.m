function my_spmbatch_run_PSCanalysis(sub,ses,run,task,datpath,paramsfile)

load(paramsfile)

global spmpath
spmpath = params.spm_path;

if params.onVSC, [datpath,params] = before_run_pscaVSC(datpath,sub,ses,task,params); end

try
    %% make batch
    params = my_spmbatch_pscanalysis(sub,ses,run,task,datpath,params);

catch e
    fprintf(['\nError processing ' num2str(sub,['%0' num2str(params.sub_digits) 'd']) ' ses-' num2str(ses,'%03d') ' run-' num2str(run,'%02d') ' task-' task '\n']);

    nlogfname = fullfile(datpath,['error_fmri_PSC_' num2str(sub,['%0' num2str(params.sub_digits) 'd']) '_ses-' num2str(ses,'%03d') '_run-' num2str(run,'%02d') '_task-' task '.txt']);

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

%----------------------------------------------------------------------
function [datpath,params] = before_run_pscaVSC(datpath,sub,ses,task,params)

fprintf('\nStart copying results\n');

orig_datpath = datpath;
pathsplit = split(datpath,'/data/');
new_datpath = ['/scratch/' pathsplit{end}];

try
    if ~isfolder(new_datpath), mkdir(new_datpath); end

    substring = ['sub-' num2str(sub,['%0' num2str(params.sub_digits) 'd'])];

    sesstring = ['ses-' num2str(ses,'%02d')];
    if ~isfolder(fullfile(datpath,substring,sesstring)), sesstring = ['ses-' num2str(ses,'%03d')]; end

    orig_subpath = fullfile(orig_datpath,substring,sesstring);
    new_subpath = fullfile(new_datpath,substring,sesstring);

    if ~isfolder(orig_subpath)
        orig_subpath = fullfile(datpath,substring); 
        new_subpath = fullfile(new_datpath,substring); 
    end

    if ~isfolder(new_subpath), mkdir(new_subpath); end

    params.resultmap = ['SPMMAT-' task '_' params.SPMMAT_analysisname];

    if isfolder(fullfile(orig_subpath,params.resultmap)), copyfile(fullfile(orig_subpath,params.resultmap),fullfile(new_subpath,params.resultmap)); end

    datpath = new_datpath;

    params.orig_subpath = orig_subpath;
    params.new_subpath = new_subpath;

    params.save_intermediate_results = false;
catch
    params.onVSC = false;
end