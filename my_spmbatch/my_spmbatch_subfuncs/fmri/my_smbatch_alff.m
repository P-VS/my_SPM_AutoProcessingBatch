function alffbrain = my_smbatch_alff(fdata,LF_band,ASamplePeriod)

[nDim1 nDimTimePoints]=size(fdata);

% Get the frequency index
sampleFreq = 1/ASamplePeriod;
NyF = sampleFreq/2;

if isinf(LF_band(2)), LF_band(2) = NyF; end

fdata = fdata - repmat(mean(fdata,2),[1,nDimTimePoints]);

fft_fdata = 2*abs(fft(fdata,[],2))/nDimTimePoints;
dimFFT = size(fft_fdata);
fft_fdata = fft_fdata(:,1:(dimFFT(2)/2) +1);
dimFFT = size(fft_fdata);

f = 0 : NyF/(dimFFT(2)-1) : NyF;
fincl = find(and(f >= LF_band(1),f <= LF_band(2)));

alffbrain = sum(fft_fdata(:,min(fincl):max(fincl)),2);