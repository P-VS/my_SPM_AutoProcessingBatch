function my_spmbatch_start_aslboldpreprocessing(sublist,nsessions,task,datpath,params)

params.func.isaslbold = true;
params.func.meepi = true;
params.echoes = params.func.echoes;
params.do_denoising = false;

if ~contains(params.func.combination,'none'), params.func.do_echocombination = true; else params.func.do_echocombination = false; end
if contains(params.asl.splitaslbold,'meica')
    params.denoise.do_mot_derivatives = true;
    params.denoise.do_aCompCor = false;
    params.denoise.Ncomponents = 5; %if in range [0 1] then the number of aCompCor components is equal to the number of components that explain the specified percentage of variation in the signal (default=5)
    params.denoise.do_bpfilter = false;
    params.denoise.bpfilter = [0.008 Inf]; %no highpass filter is first 0, no lowpass filter is last Inf, default=[0.008 Inf]
    params.denoise.polort = 1; %order of the polynomial function used to remove the signal trend (0: only mean, 1: linear trend, 2: quadratic trend, default=2)
    params.denoise.do_ICA_AROMA = true;
    params.denoise.do_noiseregression = true;
    params.denoise.do_DUNE = false;
end
if contains(params.asl.splitaslbold,'dune')
    params.denoise.do_mot_derivatives = true;
    params.denoise.do_aCompCor = false;
    params.denoise.Ncomponents = 5; %if in range [0 1] then the number of aCompCor components is equal to the number of components that explain the specified percentage of variation in the signal (default=5)
    params.denoise.do_bpfilter = false;
    params.denoise.bpfilter = [0.008 Inf]; %no highpass filter is first 0, no lowpass filter is last Inf, default=[0.008 Inf]
    params.denoise.polort = 1; %order of the polynomial function used to remove the signal trend (0: only mean, 1: linear trend, 2: quadratic trend, default=2)
    params.denoise.do_ICA_AROMA = false;
    params.denoise.do_noiseregression = false;
    params.denoise.do_DUNE = true;
    params.denoise.DUNE_part = 'bold';
end
if contains(params.asl.splitaslbold,'filter')
    params.func.denoise = true;
    params.denoise.do_mot_derivatives = true;
    params.denoise.do_aCompCor = false;
    params.denoise.Ncomponents = 5; %if in range [0 1] then the number of aCompCor components is equal to the number of components that explain the specified percentage of variation in the signal (default=5)
    params.denoise.do_bpfilter = true;
    params.denoise.bpfilter = [0.008 0.1]; %no highpass filter is first 0, no lowpass filter is last Inf, default=[0.008 Inf]
    params.denoise.polort = 1; %order of the polynomial function used to remove the signal trend (0: only mean, 1: linear trend, 2: quadratic trend, default=2)
    params.denoise.do_ICA_AROMA = false;
    params.denoise.do_noiseregression = true;
    params.denoise.do_DUNE = false;
    params.denoise.DUNE_part = 'bold';
end

if ~isfield(params,'onVSC'), params.onVSC=false; end

if params.onVSC
    params.save_intermediate_results = false;
    params.loadmaxvols = 1500;
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

                    my_spmbatch_run_aslboldpreprocessing(datlist(i,1),datlist(i,2),datlist(i,3),task{k},datpath,fullfile(datpath,paramsfile));
        
                    fprintf(['Done preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])
                end

                delete(gcp("nocreate"));
            end
        else
            for i=1:numel(datlist(:,1))
                itstart = tic;
    
                my_spmbatch_run_aslboldpreprocessing(datlist(i,1),datlist(i,2),datlist(i,3),task{k},datpath,fullfile(datpath,paramsfile));
    
                itstop = toc(itstart);
    
                fprintf(['subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
            end
        end
    end
end

delete(fullfile(datpath,paramsfile))