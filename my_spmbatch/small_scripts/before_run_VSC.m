function [datpath,params] = before_run_VSC(datpath,sub,ses,params)

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

    if isfolder(fullfile(orig_subpath,'anat')), copyfile(fullfile(orig_subpath,'anat'),fullfile(new_subpath,'anat')); end
    if isfolder(fullfile(orig_subpath,'fmap')), copyfile(fullfile(orig_subpath,'fmap'),fullfile(new_subpath,'fmap')); end
    if isfolder(fullfile(orig_subpath,'func')), copyfile(fullfile(orig_subpath,'func'),fullfile(new_subpath,'func')); end
    if isfolder(fullfile(orig_subpath,'perf')), copyfile(fullfile(orig_subpath,'perf'),fullfile(new_subpath,'perf')); end

    if isfield(params,'preprocfmridir') && isfolder(fullfile(orig_subpath,params.preprocfmridir)), copyfile(fullfile(orig_subpath,params.preprocfmridir),fullfile(new_subpath,params.preprocfmridir)); end

    datpath = new_datpath;

    params.orig_subpath = orig_subpath;
    params.new_subpath = new_subpath;
    
    params.save_intermediate_results = false;
catch
    params.onVSC = false;
end