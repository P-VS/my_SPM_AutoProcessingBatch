function my_spmbatch_start_vbmprocessing(sublist,nsessions,datpath,params)

if params.onVSC
    params.save_intermediate_results = false;
end

t = datetime('now','Format','yyMMddHHmmss');
paramsfile = ['params_' char(t) '.mat'];
save(fullfile(datpath,paramsfile),'params')

datlist = zeros(numel(sublist)*numel(nsessions),2);

dpos = 1;
for i = 1:numel(sublist)
    for j = 1:numel(nsessions)
        datlist(dpos,1) = sublist(i);
        datlist(dpos,2) = nsessions(j);

        dpos = dpos+1;
    end
end

if params.use_parallel
    for j=1:numpacks
        if (j*params.maxprocesses)<=numel(datlist(:,1))
            maxruns = params.maxprocesses;
        else
            maxruns = params.maxprocesses-((j*params.maxprocesses)-numel(datlist(:,1)));
        end

        parfor is = 1:maxruns
            i = (j-1)*params.maxprocesses+is;

            fprintf(['Start VBM for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])

            my_spmbatch_run_vbmpreprocessing(datlist(i,1),datlist(i,2),datpath,fullfile(datpath,paramsfile));
            
            fprintf(['Done VBM for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) '\n'])
        end

        delete(gcp("nocreate"));
    end
else
    for i=1:numel(datlist(:,1))
        itstart = tic;

        my_spmbatch_run_vbmpreprocessing(datlist(i,1),datlist(i,2),datpath,fullfile(datpath,paramsfile));

        itstop = toc(itstart);

        fprintf(['subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
    end
end

delete(fullfile(datpath,paramsfile))