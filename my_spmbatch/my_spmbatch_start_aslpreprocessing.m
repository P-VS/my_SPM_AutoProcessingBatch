function my_spmbatch_start_aslpreprocessing(sublist,nsessions,datpath,params)

if ~params.asl.mruns, params.asl.runs = [1]; end

t = datetime('now','Format','yyMMddHHmmss');
paramsfile = ['params_' char(t) '.mat'];
save(fullfile(datpath,paramsfile),'params')

datlist = zeros(numel(sublist)*numel(nsessions)*numel(params.asl.runs),3);

for kt = 1:numel(params.asl.runs)
    dpos = 1;
    for i = 1:numel(sublist)
        for j = 1:numel(nsessions)
                datlist(dpos,1) = sublist(i);
                datlist(dpos,2) = nsessions(j);
                datlist(dpos,3) = params.asl.runs(kt);
        
                dpos = dpos+1;
            end
    end
    
    numpacks = ceil(numel(datlist(:,1))/params.maxprocesses);
    
    for i=1:numel(datlist(:,1))
        itstart = tic;

        my_spmbatch_run_aslpreprocessing(datlist(i,1),datlist(i,2),datlist(i,3),datpath,fullfile(datpath,paramsfile));

        itstop = toc(itstart);

        fprintf(['subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' run ' num2str(datlist(i,3)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
    end
end

delete(fullfile(datpath,paramsfile))