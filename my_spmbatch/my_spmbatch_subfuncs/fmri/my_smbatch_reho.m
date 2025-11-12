function rehobrain = my_smbatch_reho(AllVolume,Vfunc,mask,NVoxel)

% Info on the approach based on the reho: 
% Zang YF, Jiang TZ, Lu YL, He Y and Tian LX, Regional Homogeneity Approach
% to fMRI Data Analysis. NeuroImage, Vol.22, No.1, 2004, 394-400.

[nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
M=nDim1;N=nDim2;O=nDim3;
isize = [nDim1 nDim2 nDim3]; 
vsize = sqrt(sum(Vfunc(1).mat(1:3,1:3).^2));

AllVolume=reshape(AllVolume,[],nDimTimePoints)';

CUTNUMBER = 10;

%% Calcualte the rank

fprintf('\n\t Rank calculating...');

Ranks_AllVolume = repmat(zeros(1,size(AllVolume,2)), [size(AllVolume,1), 1]);

SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER);
for iCut=1:CUTNUMBER
    if iCut~=CUTNUMBER
        Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
    else
        Segment = (iCut-1)*SegmentLength+1 : size(AllVolume,2);
    end
    
    AllVolume_Piece = AllVolume(:,Segment);
    nVoxels_Piece = size(AllVolume_Piece,2);
    
    [AllVolume_Piece,SortIndex] = sort(AllVolume_Piece,1);
    db=diff(AllVolume_Piece,1,1);
    db = db == 0;
    sumdb=sum(db,1);
    
    SortedRanks = repmat([1:nDimTimePoints]',[1,nVoxels_Piece]);
    % For those have the same values at the current time point and previous time point (ties)
    if any(sumdb(:))
        TieAdjustIndex=find(sumdb);
        for i=1:length(TieAdjustIndex)
            ranks=SortedRanks(:,TieAdjustIndex(i));
            ties=db(:,TieAdjustIndex(i));
            tieloc = [find(ties); nDimTimePoints+2];
            maxTies = numel(tieloc);
            tiecount = 1;
            while (tiecount < maxTies)
                tiestart = tieloc(tiecount);
                ntied = 2;
                while(tieloc(tiecount+1) == tieloc(tiecount)+1)
                    tiecount = tiecount+1;
                    ntied = ntied+1;
                end
                % Compute mean of tied ranks
                ranks(tiestart:tiestart+ntied-1) = ...
                    sum(ranks(tiestart:tiestart+ntied-1)) / ntied;
                tiecount = tiecount + 1;
            end
            SortedRanks(:,TieAdjustIndex(i))=ranks;
        end
    end
    clear db sumdb;
    SortIndexBase = repmat([0:nVoxels_Piece-1].*nDimTimePoints,[nDimTimePoints,1]);
    SortIndex=SortIndex+SortIndexBase;
    clear SortIndexBase;
    Ranks_Piece = zeros(nDimTimePoints,nVoxels_Piece);
    Ranks_Piece(SortIndex)=SortedRanks;
    clear SortIndex SortedRanks;
    
    Ranks_AllVolume(:,Segment) = Ranks_Piece;
    fprintf('.');
end

Ranks_AllVolume = reshape(Ranks_AllVolume,[nDimTimePoints,nDim1 nDim2 nDim3]);



%% calulate the kcc for the data set
% ------------------------------------------------------------------------
fprintf('\n\t Calculate the kcc on voxel by voxel for the data set.');
K = zeros(M,N,O); 
switch NVoxel
    case 27  
        for i = 2:M-1
            for j = 2:N-1
                for k = 2:O-1                            
                    block = Ranks_AllVolume(:,i-1:i+1,j-1:j+1,k-1:k+1);
                    mask_block = mask(i-1:i+1,j-1:j+1,k-1:k+1);
                    if mask_block(2,2,2)~=0
                        %YWe also calculate the ReHo value of the voxels near the border of the brain mask, users should be cautious when dicussing the results near the border.
                        R_block=reshape(block,[],27);
                        mask_R_block = R_block(:,reshape(mask_block,1,27) > 0);
                        K(i,j,k) = f_kendall(mask_R_block);  
                    end %end if	
                end%end k 
            end% end j			
			fprintf('.');
        end%end i
        fprintf('\t The reho of the data set was finished.\n');
    case 19  
        mask_cluster_19=ones(3,3,3);
        mask_cluster_19(1,1,1) = 0;    mask_cluster_19(1,3,1) = 0;    mask_cluster_19(3,1,1) = 0;    mask_cluster_19(3,3,1) = 0;
        mask_cluster_19(1,1,3) = 0;    mask_cluster_19(1,3,3) = 0;    mask_cluster_19(3,1,3) = 0;    mask_cluster_19(3,3,3) = 0;
        %The element in the mask could be 1 other than 127. Fixed the bug of computing ReHo with 7 voxels or 19 voxels.
        for i = 2:M-1
            for j = 2:N-1
                for k = 2:O-1                            
                    block = Ranks_AllVolume(:,i-1:i+1,j-1:j+1,k-1:k+1);
                    mask_block = mask(i-1:i+1,j-1:j+1,k-1:k+1);
                    if mask_block(2,2,2)~=0
                        %We also calculate the ReHo value of the voxels near the border of the brain mask, users should be cautious when dicussing the results near the border.
                        mask_block=mask_block.*mask_cluster_19;
                        %The element in the mask could be 1 other than 127. Fixed the bug of computing ReHo with 7 voxels or 19 voxels.
                        R_block=reshape(block,[],27);
                        mask_R_block = R_block(:,reshape(mask_block,1,27) > 0); %The element in the mask could be 1 other than 127. Fixed the bug of computing ReHo with 7 voxels or 19 voxels. %> 2);
                        K(i,j,k) = f_kendall(mask_R_block);  
                    end%end if
                end%end k
            end%end j	
			fprintf('.');
        end%end i
        fprintf('\t The reho of the data set was finished.\n');
    case 7   
        mask_cluster_7=ones(3,3,3);
        mask_cluster_7(1,1,1) = 0;    mask_cluster_7(1,2,1) = 0;     mask_cluster_7(1,3,1) = 0;      mask_cluster_7(1,1,2) = 0;
        mask_cluster_7(1,3,2) = 0;    mask_cluster_7(1,1,3) = 0;     mask_cluster_7(1,2,3) = 0;      mask_cluster_7(1,3,3) = 0;
        mask_cluster_7(2,1,1) = 0;    mask_cluster_7(2,3,1) = 0;     mask_cluster_7(2,1,3) = 0;      mask_cluster_7(2,3,3) = 0;
        mask_cluster_7(3,1,1) = 0;    mask_cluster_7(3,2,1) = 0;     mask_cluster_7(3,3,1) = 0;      mask_cluster_7(3,1,2) = 0;
        mask_cluster_7(3,3,2) = 0;    mask_cluster_7(3,1,3) = 0;     mask_cluster_7(3,2,3) = 0;      mask_cluster_7(3,3,3) = 0;
        %The element in the mask could be 1 other than 127. Fixed the bug of computing ReHo with 7 voxels or 19 voxels.
        for i = 2:M-1
            for j = 2:N-1
                for k = 2:O-1
                    block = Ranks_AllVolume(:,i-1:i+1,j-1:j+1,k-1:k+1);
                    mask_block = mask(i-1:i+1,j-1:j+1,k-1:k+1);
                    if mask_block(2,2,2)~=0
                        %We also calculate the ReHo value of the voxels near the border of the brain mask, users should be cautious when dicussing the results near the border.
                        mask_block=mask_block.*mask_cluster_7;
                        %The element in the mask could be 1 other than 127. Fixed the bug of computing ReHo with 7 voxels or 19 voxels.
                        R_block=reshape(block,[],27);
                        mask_R_block = R_block(:,reshape(mask_block,1,27) > 0); %The element in the mask could be 1 other than 127. Fixed the bug of computing ReHo with 7 voxels or 19 voxels. %> 2);
                        K(i,j,k) = f_kendall(mask_R_block);  
                    end%end if
                end%end k
            end%end j		
			fprintf('.');
        end%end i
        fprintf('\t The reho of the data set was finished.\n');
    otherwise
        error('The second parameter should be 7, 19 or 27. Please re-exmamin it.');
end %end switch
rehobrain = K;

clear K

end

% calculate kcc for a time series
%---------------------------------------------------------------------------
function B = f_kendall(A)
nk = size(A); n = nk(1); k = nk(2);
SR = sum(A,2); SRBAR = mean(SR);
S = sum(SR.^2) - n*SRBAR^2;
B = 12*S/k^2/(n^3-n);
end