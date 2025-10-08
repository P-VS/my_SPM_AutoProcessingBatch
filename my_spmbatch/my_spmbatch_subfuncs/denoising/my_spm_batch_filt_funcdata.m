function funcdat = my_spm_batch_filt_funcdata(funcdat, ltrh, htrh, tr)

Fs = 1/tr; %Sample frequency
Ny = Fs/2; %Nyquist frequency

if ltrh==0 || ltrh<0.0, ltrh=-1; end
if htrh==Ny || htrh>Ny, htrh=-1; end

s = size(funcdat);

mean_funcdat = mean(funcdat,2);
funcdat = funcdat - repmat(mean_funcdat,[1,s(end)]);

nvox = numel(funcdat(:,1));
for i=1:5000:nvox
    iend = min(i+5000-1,nvox);

    if ~(ltrh==-1)
        bpdat = brant_Filter_FFT_Butterworth(funcdat(i:iend,:)', -1, ltrh, Fs, 0)';
        funcdat(i:iend,:) = funcdat(i:iend,:)-bpdat;
    end
    if ~(htrh==-1), funcdat(i:iend,:) = brant_Filter_FFT_Butterworth(funcdat(i:iend,:)', -1, htrh, Fs, 0)'; end
end

funcdat = funcdat + repmat(mean_funcdat,[1,s(end)]);

clear mean_funcdat