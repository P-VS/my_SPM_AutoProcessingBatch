function params = my_spmbatch_fmrilevel1processing(sub,ses,run,task,datpath,params)

%% Search for the data folders

[ppparams,params,datpath] = my_spmbatch_1stlevel_FindData(sub,ses,run,task,datpath,params);

%% fMRI model specification

if params.func.mruns && contains(params.func.use_runs,'separately')
    ppparams.resultfolder = ['SPMMAT-' task '_' params.analysisname '_run-' num2str(run)];
    ppparams.resultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.analysisname '_run-' num2str(run)]);
else
    ppparams.resultfolder = ['SPMMAT-' task '_' params.analysisname];
    ppparams.resultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.analysisname]);
end

params.resultmap = ppparams.resultfolder;

if exist(ppparams.resultmap,'dir'); rmdir(ppparams.resultmap,'s'); end
mkdir(ppparams.resultmap)

[fmri_spec,ppparams] = my_spmbatch_1stlevel_DefineModel(sub,ses,run,task,datpath,params,ppparams);

spm_run_fmri_spec(fmri_spec);

SPM_file = fullfile(ppparams.resultmap,'SPM.mat');

%% Optimize GLM with TEDM

if params.optimize_HRF && contains(params.analysis_type,'GLM')
    load(SPM_file)
    Sess = SPM.Sess;

    SPM_file = my_spmmbatch_tedm(SPM_file,ppparams.resultmap,ppparams.mask_file); 

    load(SPM_file)
    SPM.Sess = Sess;
    save(SPM_file,'SPM')
end

%% Model estimation

fmri_est.spmmat(1) = {SPM_file};
fmri_est.write_residuals = 0;
fmri_est.method.Classical = 1;
  
spm_run_fmri_est(fmri_est);

%% Contrast Manager

con.spmmat(1) = {SPM_file};

for ic=1:numel(params.contrast)
    contrastname='';
    weights = [];

    for ir=1:numel(params.iruns)
        switch params.analysis_type
            case 'GLM'
                if ~params.add_derivatives; ncondcol = 1; else ncondcol = 3; end
                if params.add_parametricModulation, ncondcol=ncondcol+numel(ppparams.edat{ir}.weight); end
            case 'FIR'
                ncondcol = fmri_spec.bases.fir.order;
        end

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
                if strcmp(lower(params.contrast(ic).conditions{icn}),lower(ppparams.edat{ir}.conditions{icn2}.name)) 
                    switch params.analysis_type
                        case 'GLM'
                            indx=(icn2-1)*ncondcol+1; 
                        case 'FIR'
                            indx=(icn2-1)*ncondcol+1:icn2*ncondcol; 
                    end
                end
            end
    
            if indx(1)>0; subweights(indx)=params.contrast(ic).vector(icn); end
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

    con.consess{ic}.tcon.name = contrastname;
    con.consess{ic}.tcon.weights = weights;
    
    con.consess{ic}.tcon.sessrep = 'none';
end

con.delete = 0;

spm_run_con(con);

%% SPM Results

if params.save_spm_results
    results.spmmat = {SPM_file};
    results.conspec.titlestr = '';
    results.conspec.contrasts = Inf;
    results.conspec.threshdesc = params.threshold_correction;
    results.conspec.thresh = params.pthreshold;
    results.conspec.extent = params.kthreshold;
    results.conspec.conjunction = 1;
    results.conspec.mask.none = 1;
    results.units = 1;

    oi = 1;
    if params.save_thresholded_map
        results.export{oi}.tspm.basename = 'thres';
        oi = oi+1;
    end
    if params.save_binary_mask
        results.export{oi}.binary.basename = 'bin';
        oi = oi+1;
    end
    if params.save_naray
        results.export{oi}.nary.basename = 'n-aray';
        oi = oi+1;
    end
    if params.save_csv_file
        results.export{oi}.csv = true;
        oi = oi+1;
    end
    if params.save_pdf_file
        results.export{oi}.pdf = true;
        oi = oi+1;
    end
    if params.save_tiff_file
        results.export{oi}.tif = true;
        oi = oi+1;
    end

    my_spmbatch_run_results(results);
end
