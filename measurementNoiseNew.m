function [R, varargout] = measurementNoiseNew(xseg,fs)
%new method of calculating measurement noise variance based on PSD

numFrame = size(xseg,1);
noise_cov = zeros(1,numFrame);
spectral_flatness = zeros(1,numFrame);
silent_inds = zeros(1,numFrame);

%order estimation for voiced and silent frames
for k = 1:numFrame
    
    [c, lag] = xcorr(xseg(k,:),'coeff');
    %calculating power spectral density from ACF
    psd = (fftshift(abs(fft(c))));
    psd = psd(round(length(psd)/2):end);
    freq = (fs * (0:length(c)/2))/length(c);
    %keeping positive lags only since ACF is symmetrical
    c = c(find(lag == 0):length(c));
    lag = lag(find(lag == 0):length(lag));
    %keep frequencies from 100Hz to 2kHz
    freq_2kHz = find(freq>= 100 & freq<=2000);
    psd_2kHz = psd(freq_2kHz);
    spectral_flatness(k) = geomean(psd_2kHz)/mean(psd_2kHz);
    
end

normalized_flatness = spectral_flatness/max(spectral_flatness);
threshold = 0.707;
for k = 1:numFrame
    if normalized_flatness(k) >= threshold
        noise_cov(k) = var(xseg(k,:));
        silent_inds(k) = 1;
    end
end
R = max(noise_cov);
varargout{1} = silent_inds;
end

