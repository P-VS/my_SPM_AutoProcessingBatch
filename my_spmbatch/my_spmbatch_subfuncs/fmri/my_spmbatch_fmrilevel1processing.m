function matlabbatch = my_spmbatch_fmrilevel1processing(sub,ses,run,task,datpath,params)

matlabbatch = {};

%% Search for the data folders
ppparams = my_spmbatch_1stlevel_FindData(sub,ses,run,task,datpath,params);

%% fMRI model specification
if params.func.mruns && contains(params.func.use_runs,'separately')
    ppparams.resultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.analysisname '_run-' num2str(run)]);
else
    ppparams.resultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.analysisname]);
end

if exist(ppparams.resultmap,'dir'); rmdir(ppparams.resultmap,'s'); end
mkdir(ppparams.resultmap)

[matlabbatch,ppparams] = my_spmbatch_1stlevel_DefineModel(sub,ses,run,task,datpath,params,ppparams,matlabbatch);

%% Optimize GLM with TEDM
if params.optimize_HRF
    spm_jobman('run', matlabbatch)

    clear matlabbatch

    SPM_file = fullfile(ppparams.resultmap,'SPM.mat');
    SPM_file = my_spmmbatch_tedm(SPM_file,ppparams.resultmap,ppparams.mask_file);

    mbidx = 1;
    matlabbatch{mbidx}.spm.stats.fmri_est.spmmat(1) = {SPM_file};
else
    mbidx = 2;
    matlabbatch{mbidx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{mbidx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
end

%% Model estimation

matlabbatch{mbidx}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{mbidx}.spm.stats.fmri_est.method.Classical = 1;
  
mbidx=mbidx+1;

%% Contrast Manager

matlabbatch{mbidx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{mbidx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

for ic=1:numel(params.contrast)
    contrastname='';
    weights = [];

    for ir=1:numel(params.iruns)
        if ~params.add_derivatives; ncondcol = 1; else ncondcol = 3; end
        if params.add_parametricModulation, ncondcol=ncondcol+numel(ppparams.edat{ir}.weight); end

        subweights = zeros(1,ncondcol*numel(numel(ppparams.edat{ir}.onset)));

        if params.add_regressors
            rpdat = load(ppparams.frun(ir).confoundsfile);
            subweights=[subweights zeros(1,numel(rpdat(1,:)))];
        end
        for icn=1:numel(params.contrast(ic).conditions)
            if params.contrast(ic).vector(icn)>0; contrastname = [contrastname ' + ' params.contrast(ic).conditions{icn}]; end
            if params.contrast(ic).vector(icn)<0; contrastname = [contrastname ' - ' params.contrast(ic).conditions{icn}]; end
        
            indx=0;
            for icn2=1:numel(ppparams.edat{ir}.conditions)
                if strcmp(lower(params.contrast(ic).conditions{icn}),lower(ppparams.edat{ir}.conditions{icn2}.name)); indx=(icn2-1)*ncondcol+1; end
            end
    
            if indx>0; subweights(indx)=params.contrast(ic).vector(icn); end
        end

        subweights = repmat(subweights,1,numel(params.func.echoes));

        weights=[weights subweights];
    end

    postmp = find(weights>0);
    if ~isempty(postmp)
        weights(postmp) = weights(postmp)/numel(postmp);
    end
    negtmp = find(weights<0);
    if ~isempty(negtmp)
        weights(negtmp) = weights(negtmp)/numel(negtmp);
    end

    matlabbatch{mbidx}.spm.stats.con.consess{ic}.tcon.name = contrastname;
    matlabbatch{mbidx}.spm.stats.con.consess{ic}.tcon.weights = weights;
    
    matlabbatch{mbidx}.spm.stats.con.consess{ic}.tcon.sessrep = 'none';
end

matlabbatch{mbidx}.spm.stats.con.delete = 0;

mbidx=mbidx+1;

%% SPM Results
if params.save_spm_results
    matlabbatch{mbidx}.spm.stats.results.spmmat = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{mbidx-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{mbidx}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{mbidx}.spm.stats.results.conspec.contrasts = Inf;
    matlabbatch{mbidx}.spm.stats.results.conspec.threshdesc = params.threshold_correction;
    matlabbatch{mbidx}.spm.stats.results.conspec.thresh = params.pthreshold;
    matlabbatch{mbidx}.spm.stats.results.conspec.extent = params.kthreshold;
    matlabbatch{mbidx}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{mbidx}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{mbidx}.spm.stats.results.units = 1;

    oi = 1;
    if params.save_thresholded_map
        matlabbatch{mbidx}.spm.stats.results.export{oi}.tspm.basename = 'thres';
        oi = oi+1;
    end
    if params.save_binary_mask
        matlabbatch{mbidx}.spm.stats.results.export{oi}.binary.basename = 'bin';
        oi = oi+1;
    end
    if params.save_naray
        matlabbatch{mbidx}.spm.stats.results.export{oi}.nary.basename = 'n-aray';
        oi = oi+1;
    end
    if params.save_csv_file
        matlabbatch{mbidx}.spm.stats.results.export{oi}.csv = true;
        oi = oi+1;
    end
    if params.save_pdf_file
        matlabbatch{mbidx}.spm.stats.results.export{oi}.pdf = true;
        oi = oi+1;
    end
    if params.save_tiff_file
        matlabbatch{mbidx}.spm.stats.results.export{oi}.tif = true;
        oi = oi+1;
    end
end
