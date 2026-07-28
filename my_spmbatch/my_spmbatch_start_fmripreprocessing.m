function my_spmbatch_start_fmripreprocessing(sublist,nsessions,task,datpath,params)

params.func.isaslbold = false;
params.preprocess_asl = false;
params.asl.splitaslbold = 'none';

if ~params.func.mruns, params.func.runs = [1]; end
if params.func.meepi && numel(params.func.echoes)==1, params.func.combination='none'; end
if params.func.meepi && ~contains(params.func.combination,'none'), params.func.do_echocombination = true; else params.func.do_echocombination = false; end

if params.denoise.do_DUNE
    params.denoise.do_noiseregression = false;
    params.denoise.do_ICA_AROMA = false;
    params.denoise.do_aCompCor = true;
    params.denoise.do_bpfilter = true;
    params.denoise.bpfilter = [0.008 Inf];
    params.denoise.do_mot_derivatives = true;
end

if params.denoise.do_ICA_AROMA
    params.denoise.do_noiseregression = true;
    params.denoise.do_DUNE = false;
    params.denoise.do_aCompCor = true;
    params.denoise.do_bpfilter = true;
    params.denoise.bpfilter = [0.008 Inf];
    params.denoise.do_mot_derivatives = true;
end

if params.denoise.do_aCompCor
    params.denoise.do_noiseregression = true;
    params.denoise.do_bpfilter = true;
    params.denoise.bpfilter = [0.008 Inf];
    params.denoise.do_mot_derivatives = true;
end

if ~isfield(params,'onVSC'), params.onVSC=false; end

if params.onVSC
    params.save_intermediate_results = false;
    params.loadmaxvols = 800;
end

t = datetime('now','Format','yyMMddHHmmss');
paramsfile = ['params_' char(t) '.mat'];
save(fullfile(datpath,paramsfile),'params')

datlist = zeros(numel(sublist)*numel(nsessions),3);

for kt = 1:numel(params.func.runs)
    for k = 1:numel(task)
        dpos = 1;
        for i = 1:numel(sublist)
            for j = 1:numel(nsessions)
                datlist(dpos,1) = sublist(i);
                datlist(dpos,2) = nsessions(j);
                datlist(dpos,3) = params.func.runs(kt);
        
                dpos = dpos+1;
            end
        end
        
        numpacks = ceil(numel(datlist(:,1))/params.maxprocesses);
        
        if params.use_parallel
            for j=1:numpacks
                if (j*params.maxprocesses)<=numel(datlist(:,1))
                    maxruns = params.maxprocesses;
                else
                    maxruns = params.maxprocesses-((j*params.maxprocesses)-numel(datlist(:,1)));
                end

                parfor is = 1:maxruns
                    i = (j-1)*params.maxprocesses+is;

                    fprintf(['Start preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])
        
                    my_spmbatch_run_fmripreprocessing(datlist(i,1),datlist(i,2),datlist(i,3),task{k},datpath,fullfile(datpath,paramsfile));
        
                    fprintf(['Done preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])
                end

                delete(gcp("nocreate"));
            end
        else
            for i=1:numel(datlist(:,1))
                itstart = tic;
    
                my_spmbatch_run_fmripreprocessing(datlist(i,1),datlist(i,2),datlist(i,3),task{k},datpath,fullfile(datpath,paramsfile));
    
                itstop = toc(itstart);
    
                fprintf(['subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
            end
        end
    end
end

delete(fullfile(datpath,paramsfile))