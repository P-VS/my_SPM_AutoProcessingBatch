function my_spmbatch_start_PSCanalysis(sublist,nsessions,datpath,params)

if params.onVSC || params.use_parallel, params.save_spm_results = false; end

if ~isfield(params,'onVSC'), params.onVSC=false; end

if params.onVSC
    params.use_parallel = false;
    params.save_intermediate_results = false;
    params.loadmaxvols = 1000;
end

t = datetime('now','Format','yyMMddHHmmss');
paramsfile = ['params_' char(t) '.mat'];
save(fullfile(datpath,paramsfile),'params')

datlist = zeros(numel(sublist)*numel(nsessions),3);

dpos = 1;
for i = 1:numel(sublist)
    for j = 1:numel(nsessions)
        for k = 1:numel(params.func.runs)
            datlist(dpos,1) = sublist(i);
            datlist(dpos,2) = nsessions(j);
            datlist(dpos,3) = params.func.runs(k);
    
            dpos = dpos+1;
        end
    end
end

numpacks = ceil(numel(datlist(:,1))/params.maxprocesses);

for k = 1:numel(params.task)
    if params.use_parallel
        for j=1:numpacks
            if (j*params.maxprocesses)<=numel(datlist(:,1))
                maxruns = params.maxprocesses;
            else
                maxruns = params.maxprocesses-((j*params.maxprocesses)-numel(datlist(:,1)));
            end

            parfor is = 1:maxruns
                i = (j-1)*params.maxprocesses+is;

                fprintf(['Start PSC for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])

                my_spmbatch_run_PSCanalysis(datlist(i,1),datlist(i,2),datlist(i,3),params.task{k},datpath,fullfile(datpath,paramsfile));
    
                fprintf(['Done PSC for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])
            end
        end
    else  
        for i=1:numel(datlist(:,1))
            fprintf(['Start PSC for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])

            my_spmbatch_run_PSCanalysis(datlist(i,1),datlist(i,2),datlist(i,3),params.task{k},datpath,fullfile(datpath,paramsfile));

            fprintf(['Done PSC for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])
        end
    end
end

delete(fullfile(datpath,paramsfile))