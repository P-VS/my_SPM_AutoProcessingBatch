function ppparams = my_spmbatch_1stlevel_FindData(sub,ses,run,task,datpath,params)

%% Search for the data folders

ppparams.substring = ['sub-' num2str(sub,['%0' num2str(params.sub_digits) 'd'])];

ppparams.sesstring = ['ses-' num2str(ses,'%02d')];
if ~isfolder(fullfile(datpath,ppparams.substring,ppparams.sesstring)), ppparams.sesstring = ['ses-' num2str(ses,'%03d')]; end

ppparams.subpath = fullfile(datpath,ppparams.substring,ppparams.sesstring);

if ~isfolder(ppparams.subpath), ppparams.subpath = fullfile(datpath,ppparams.substring); end

if ~isfolder(ppparams.subpath)
    fprintf(['No data folder for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.subfuncdir = fullfile(ppparams.subpath,'func');

if ~isfolder(ppparams.subfuncdir)
    fprintf(['No func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.preprocfmridir = fullfile(ppparams.subpath,params.preprocfmridir);

if ~isfolder(ppparams.subfuncdir)
    fprintf(['No preprocessed func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

for ir=1:numel(params.iruns)
    %% Search for the data files
    
    namefilters(1).name = ppparams.substring;
    namefilters(1).required = true;
    
    namefilters(2).name = ppparams.sesstring;
    namefilters(2).required = false;
    
    switch params.func.use_runs
        case 'separately'
            namefilters(3).name = ['run-' num2str(run)];
        case 'together'
            namefilters(3).name = ['run-' num2str(params.iruns(ir))];
    end
    if params.func.mruns, namefilters(3).required = true; else namefilters(3).required = false; end
    
    namefilters(4).name = ['task-' task];
    namefilters(4).required = true;
    
    %fMRI data
    
    if contains(params.modality,'fmri'), namefilters(5).name = '_bold'; end
    if contains(params.modality,'fasl') && contains(params.whichfile,'cbf'), namefilters(5).name = '_cbf'; end
    if contains(params.modality,'fasl') && contains(params.whichfile,'asl'), namefilters(5).name = '_asl'; end
    namefilters(5).required = true;

    namefilters(6).name = params.fmri_prefix;
    namefilters(6).required = true;

    funcniilist = my_spmbatch_dirfilelist(ppparams.preprocfmridir,'nii',namefilters,false);
    
    if isempty(funcniilist)
        fprintf(['No nii files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' params.task '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    for ie=1:numel(params.func.echoes)
        if params.func.meepi && params.use_echoes_as_sessions %Filter list based on echo number
            tmp = find(or(contains({funcniilist.name},['_echo-' num2str(params.func.echoes(ie))]),contains({funcniilist.name},['_e' num2str(params.func.echoes(ie))])));
            if isempty(tmp), edirniilist = funcniilist; else edirniilist = funcniilist(tmp); end
        else
            edirniilist = funcniilist;
        end
    
        prefixlist = split({edirniilist.name},'sub-');
        if numel(edirniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end
    
        tmp = find(strcmp(prefixlist,params.fmri_prefix));
        if ~isempty(tmp), ppparams.frun(ir).func(ie).funcfile = edirniilist(tmp).name; end
    
        if ~isfield(ppparams.frun(ir).func(ie),'funcfile')
            fprintf(['no preprocessed fmri data found for run ' num2str(ir) ' for echo ' num2str(params.func.echoes(ie)) '\n'])
            fprintf('\nPP_Error\n');
            return
        end
    
        Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(ir).func(ie).funcfile));
    
        for i=1:numel(Vfunc)
            ppparams.ppfmridat{ir}.sess{ie}.func{i,1} = [Vfunc(i).fname ',' num2str(i)];
        end
    end
    
    %confound file
    
    if params.add_regressors
        cnamefilters = namefilters;

        namefilters(5).name = '_bold';
        
        cnamefilters(6).name = params.confounds_prefix;
        cnamefilters(6).required = true;
    
        funcconlist = my_spmbatch_dirfilelist(ppparams.preprocfmridir,'txt',cnamefilters,false);
        
        if isempty(funcconlist)
            cnamefilters(5).name = '_asl';
            cnamefilters(5).required = true;
    
            funcconlist = my_spmbatch_dirfilelist(ppparams.preprocfmridir,'txt',cnamefilters,false);
        end

        if isempty(funcconlist)
            fprintf(['No confound file found for ' ppparams.substring ' ' ppparams.sesstring ' task-' task '\n'])
            fprintf('\nPP_Error\n');
            return
        end
        
        prefixlist = split({funcconlist.name},'sub-');
        if numel(funcconlist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end
        
        tmp = find(strcmp(prefixlist,params.confounds_prefix));
        if ~isempty(tmp), ppparams.frun(ir).confoundsfile = fullfile(funcconlist(tmp).folder,funcconlist(tmp).name); else ppparams.frun(ir).confoundsfile = ''; end        
    else
        ppparams.frun(ir).confoundsfile = '';
    end

    %events.tsv file
    
    enamefilters(1:4) = namefilters(1:4);
    
    enamefilters(5).name = '_events';
    enamefilters(5).required = true;
    
    functsvlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'tsv',enamefilters,false);
    
    if isempty(functsvlist)
        fprintf(['No events.tsv files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' task '\n'])
        ppparams.frun(ir).functsvfile = '';
    else
        ppparams.frun(ir).functsvfile = fullfile(functsvlist(1).folder,functsvlist(1).name);
    end
    
    %json file

    jnamefilters(1:4) = namefilters(1:4);
    
    jnamefilters(5).name = '_bold';
    jnamefilters(5).required = true;
    
    funcjsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',jnamefilters,false);
        
    if isempty(funcjsonlist)
        jnamefilters(5).name = '_asl';
        jnamefilters(5).required = true;
        
        funcjsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',jnamefilters,false);
    end

    if params.func.meepi
        jstmp = find(or(contains({funcjsonlist.name},'echo-1'),contains({funcjsonlist.name},'_e1')));
        funcjsonlist = funcjsonlist(jstmp);
    end
    
    if isempty(funcjsonlist)
        fprintf(['No json files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' task '\n'])
        fprintf('\nPP_Error\n');
        return
    end
    
    ppparams.frun(ir).funcjsonfile = fullfile(funcjsonlist(1).folder,funcjsonlist(1).name);
end
