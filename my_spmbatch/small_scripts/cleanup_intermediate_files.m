function cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params,subdir,save_folder)

if isfolder(fullfile(datpath,['sub-' num2str(sub)])), ppparams.substring = ['sub-' num2str(sub)]; end
if isfolder(fullfile(datpath,['sub-0' num2str(sub)])), ppparams.substring = ['sub-0' num2str(sub)]; end
if isfolder(fullfile(datpath,['sub-00' num2str(sub)])), ppparams.substring = ['sub-00' num2str(sub)]; end
if isfolder(fullfile(datpath,['sub-000' num2str(sub)])), ppparams.substring = ['sub-000' num2str(sub)]; end

if isfolder(fullfile(datpath,ppparams.substring,['ses-' num2str(ses)])), ppparams.sesstring = ['ses-' num2str(ses)]; end
if isfolder(fullfile(datpath,ppparams.substring,['ses-0' num2str(ses)])), ppparams.sesstring = ['ses-0' num2str(ses)]; end
if isfolder(fullfile(datpath,ppparams.substring,['ses-00' num2str(ses)])), ppparams.sesstring = ['ses-00' num2str(ses)]; end

subpath = fullfile(datpath,ppparams.substring,ppparams.sesstring);

subolddir = fullfile(subpath,subdir);
subnewdir = fullfile(subpath,save_folder);

if ~exist(subnewdir, 'dir')
    mkdir(subnewdir);
end

%% Delete intermediate files
if ~params.save_intermediate_results
    for i=1:numel(delfiles)
        if isfile(delfiles{i})
            if isfile(delfiles{i}{1}); delete(delfiles{i}{1}); end
        elseif isfolder(delfiles{i})
            rmdir(delfiles{i}{1},'s');
        end
    end
end

%% Move results
for i=1:numel(keepfiles)
    if isfile(keepfiles{i}{1})
        [opath,~,~] = fileparts(keepfiles{i}{1});
        if ~strcmp(opath,subnewdir)
            if params.save_intermediate_results, copyfile(keepfiles{i}{1},subnewdir); 
            else movefile(keepfiles{i}{1},subnewdir); end
        end
    end
end

end