function [ppparams,delfiles,keepfiles] = my_spmbacth_faslsubtraction(ppparams,params,delfiles,keepfiles)

jsondat = fileread(ppparams.func(1).jsonfile);
jsondat = jsondecode(jsondat);
tr = jsondat.RepetitionTime;

Vasl=spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix ppparams.perf(1).aslfile]));
fasldata = spm_read_vols(Vasl);

voldim = size(fasldata);

mask = my_spmbatch_mask(fasldata);

conidx = 2:2:voldim(4);
labidx = 1:2:voldim(4);

ppparams.asl.conidx = conidx;
ppparams.asl.labidx = labidx;

deltamdata = zeros([voldim(1)*voldim(2)*voldim(3),voldim(4)]);
mdeltamdata = zeros([voldim(1)*voldim(2)*voldim(3),1]);

fasldata = reshape(fasldata,[voldim(1)*voldim(2)*voldim(3),voldim(4)]);

ncondat = fasldata(mask>0,:);
nlabdat = fasldata(mask>0,:);

%Do subtraction
switch params.asl.interp
    case 'sinc'
        for p=1:voldim(4)
            if sum(conidx==p)==0
                % 6 point sinc interpolation
                idx=p+[-5 -3 -1 1 3 5];
                idx(find(idx<min(conidx)))=min(conidx);
                idx(find(idx>max(conidx)))=max(conidx);
                normloc=3.5;
            
                ncondat(:,p)=sinc_interpVec(fasldata(mask>0, idx),normloc);
            end
            if sum(labidx==p)==0
                % 6 point sinc interpolation
                idx=p+[-5 -3 -1 1 3 5];
                idx(find(idx<min(labidx)))=min(labidx);
                idx(find(idx>max(labidx)))=max(labidx);
                normloc=3.5;
  
                nlabdat(:,p)=sinc_interpVec(fasldata(mask>0, idx),normloc);
            end
        end
    case 'spline'
        ncondat=interp1((conidx-1)*tr,fasldata(mask>0, conidx)',[0:1:voldim(4)-1]*tr,'spline','extrap')';
        nlabdat=interp1((labidx-1)*tr,fasldata(mask>0, labidx)',[0:1:voldim(4)-1]*tr,'spline','extrap')';
end

if contains(params.asl.Tres,'orig')
    deltamdata = zeros([voldim(1)*voldim(2)*voldim(3),voldim(4)]);
    deltamdata(mask>0,:) = (ncondat-nlabdat);
    deltamdata(isnan(deltamdata)) = 0;
    deltamdata = deltamdata(:,2:end-1);
    voldim(4) = voldim(4)-2;
elseif contains(params.asl.Tres,'half')
    voldim(4) = numel(labidx);
    deltamdata = zeros([voldim(1)*voldim(2)*voldim(3),numel(labidx)]);
    deltamdata(mask>0,:) = (ncondat(:,labidx)-nlabdat(:,labidx));
    deltamdata(isnan(deltamdata)) = 0;
end

tmp = deltamdata(mask>0,:);
outtest = isoutlier(tmp,"median");
if ~isempty(outtest), tmp(outtest)=0; end
deltamdata(mask>0,:) = tmp;

clear tmp outtest

mdeltamdata(mask>0,:) = mean(deltamdata(mask>0,:),2);

deltamdata = reshape(deltamdata,[voldim(1),voldim(2),voldim(3),voldim(4)]);
mdeltamdata = reshape(mdeltamdata,[voldim(1),voldim(2),voldim(3)]);

clear ncondat nlabdat fasldata

nfname = split(ppparams.perf(1).aslfile,'_asl');

Vout = Vasl(1);
rmfield(Vout,'pinfo');
for iv=1:voldim(4)
    Vout.fname = fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix nfname{1} '_deltam.nii']);
    Vout.descrip = 'my_spmbatch - deltam';
    Vout.dt = [spm_type('float32'),spm_platform('bigend')];
    Vout.n = [iv 1];
    Vout = spm_write_vol(Vout,deltamdata(:,:,:,iv));
end

Vout = Vasl(1);
rmfield(Vout,'pinfo');
Vout.fname = fullfile(ppparams.subperfdir,['mean_' ppparams.perf(1).aslprefix nfname{1} '_deltam.nii']);
Vout.descrip = 'my_spmbatch - deltam';
Vout.dt = [spm_type('float32'),spm_platform('bigend')];
Vout.n = [1 1];
Vout = spm_write_vol(Vout,mdeltamdata);

clear Vout deltamdata mdeltamdata

ppparams.perf(1).meandeltam = ['mean_' ppparams.perf(1).aslprefix nfname{1} '_deltam.nii'];
ppparams.perf(1).deltamprefix = ppparams.perf(1).aslprefix;
ppparams.perf(1).deltamfile = [nfname{1} '_deltam.nii'];

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix nfname{1} '_deltam.nii'])};
delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.perf(1).meandeltam)};