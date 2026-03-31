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

if isfolder(fullfile(params.new_subpath,'perf'))
    if ~isfolder(fullfile(params.orig_subpath,'perf')), mkdir(fullfile(params.orig_subpath,'perf')); end

    dirlist = dir(fullfile(params.new_subpath,'perf')); %Make list of alll files
    
    tmp = find(strlength({dirlist.name})>4); %Remove '.' and '..'
    if ~isempty(tmp)
        dirlist = dirlist(tmp);
    
        tmp = find(~contains({dirlist.name},'._')); %Remove the hiden files from Mac from the list
        if ~isempty(tmp)
            dirlist = dirlist(tmp);
    
            tmp = find(contains({dirlist.name},'_label'));
            if ~isempty(tmp)
                for i=1:numel(tmp), copyfile(fullfile(params.new_subpath,'perf',dirlist(tmp(i)).name),fullfile(params.orig_subpath,'perf',dirlist(tmp(i)).name)); end
            end
            tmp = find(and(contains({dirlist.name},'_m0scan'),~startsWith({dirlist.name},'w')));
            if ~isempty(tmp)
                for i=1:numel(tmp), copyfile(fullfile(params.new_subpath,'perf',dirlist(tmp(i)).name),fullfile(params.orig_subpath,'perf',dirlist(tmp(i)).name)); end
            end
            tmp = find(contains({dirlist.name},'_asl'));
            if ~isempty(tmp)
                for i=1:numel(tmp), copyfile(fullfile(params.new_subpath,'perf',dirlist(tmp(i)).name),fullfile(params.orig_subpath,'perf',dirlist(tmp(i)).name)); end
            end
            tmp = find(startsWith({dirlist.name},'p0e'));
            if ~isempty(tmp)
                for i=1:numel(tmp), copyfile(fullfile(params.new_subpath,'perf',dirlist(tmp(i)).name),fullfile(params.orig_subpath,'perf',dirlist(tmp(i)).name)); end
            end
            tmp = find(startsWith({dirlist.name},'p1e'));
            if ~isempty(tmp)
                for i=1:numel(tmp), copyfile(fullfile(params.new_subpath,'perf',dirlist(tmp(i)).name),fullfile(params.orig_subpath,'perf',dirlist(tmp(i)).name)); end
            end
            tmp = find(startsWith({dirlist.name},'p2e'));
            if ~isempty(tmp)
                for i=1:numel(tmp), copyfile(fullfile(params.new_subpath,'perf',dirlist(tmp(i)).name),fullfile(params.orig_subpath,'perf',dirlist(tmp(i)).name)); end
            end
        end
    end
end

if isfolder(fullfile(params.new_subpath,'anat')), rmdir(fullfile(params.new_subpath,'anat'),'s'); end
if isfolder(fullfile(params.new_subpath,'fmap')), rmdir(fullfile(params.new_subpath,'fmap'),'s'); end
if isfolder(fullfile(params.new_subpath,'func')), rmdir(fullfile(params.new_subpath,'func'),'s'); end
if isfolder(fullfile(params.new_subpath,'perf')), rmdir(fullfile(params.new_subpath,'perf'),'s'); end

if isfield(params,'preprocfmridir')
    if isfolder(fullfile(params.new_subpath,params.preprocfmridir)), rmdir(fullfile(params.new_subpath,params.preprocfmridir),'s'); end
end
if isfield(params,'resultmap') 
    if isfolder(fullfile(params.new_subpath,params.resultmap))
        copyfile(fullfile(params.new_subpath,params.resultmap),fullfile(params.orig_subpath,params.resultmap)); 
        rmdir(fullfile(params.new_subpath,params.resultmap),'s');
    end
end

datpath = params.orig_subpath;