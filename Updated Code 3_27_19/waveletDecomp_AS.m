function [phase,pow,filtsig] = waveletDecomp_AS(freqs,signal,srate,numCycles)

%  [phase,pow,filtsig] = waveletDecomp_AS(freqs,signal,srate,numCycles)
%
%   Input:  freqs = frequencies to analyze
%           signal = signal to analyze (vector or matrix)
%           srate = sampling rate (hz)
%           numCycles = number of Morlet wavelet cycles
%   Output: phase
%           power (scaled to amplitude of input)
%           filtsig = original signal filtered at each frequency

% Check dimensions of signal
if size(signal,1) > 1
    signal = signal';
end

signal = single(signal);

% Preallocate variables for output
pow = zeros(numel(freqs),numel(signal),'single');
phase = zeros(numel(freqs),numel(signal),'single');
filtsig = zeros(numel(freqs),numel(signal),'single');

% Time for wavelet
waveTime = -2:1/srate:2;

% Length of convolution (used in fft)
% Note: The fft runs faster if the length of the transform is a power of two.
Lconv = length(waveTime) + length(signal) - 1;
Lconv2 = pow2(nextpow2(Lconv));
startIndex = ceil(length(waveTime)/2);
endIndex = Lconv - floor(length(waveTime)/2);

% Signal fft
signalFFT = fft(signal,Lconv2,2);

for i = 1:length(freqs)
    
    % Create wavelet
    % In order to create a wavelet, you need to window a sine wave with a
    % Gaussian:
    % morelet wavelet = sine wave .* gaussian window
    % sine wave = exp(2*pi*1i*f.*t)
    % gaussian window = exp(-(t^2)/(2s^2))
    % t = wave time
    % s = width of Gaussian = n/(2pi*f)
    % n = number of wavelet cycles
    s = numCycles/(2*pi*freqs(i));
    gaussianWindow = exp(-(waveTime.^2)/(2*s^2));
    sineWave = exp(2*pi*1i*freqs(i) .* waveTime);
    morletWavelet = sineWave .* gaussianWindow;
    
%     %%%%%%%%%%%%% FFT that gives the same as the convolution %%%%%%%%%%%%%%
%     % Note: The creation of the gaussian depends on a frequency; for a
%     % simple convolution with a gaussian, this frequency is not the same as
%     % the frequency you wish to filter the signal. For example A gaussian f
%     % of 30 filters your signal around 7 Hz.
%     
%     % Wavelet fft
%     waveletFFT = fft(morletWavelet,Lconv2);
% 
%     % Inverse wavelet fft
%     signalIFT = ifft(waveletFFT .* signalFFT,Lconv2);
%     signalIFT = signalIFT(1:Lconv); % remove the zero-padded signal
%     
%     % Remove the appropriate number of timepoints from the beginning and
%     % end of the time series due to convolution
%     signalIFT = signalIFT(startIndex:endIndex);
%     
%     s = numCycles/(2*pi*30);
%     gaussianWindow = (1/30) * exp(-(waveTime.^2)/(2*s^2)); % The 1/30 modifies the amplitude of the gaussian.
%     signalConv(i,:) = conv(signal,gaussianWindow,'same');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %%%%%%%%%%%%%%% Alter amplitude of Morlet Wavelet %%%%%%%%%%%%%%%%%%%%%
%     % Note: This provides a filtered signal with way to large of an
%     % amplitude; for some reason, this does not do the appropriate scaling.
%     % I got this code from morlet.m in Lindsay Vass's old analyses.
%     A = 1/sqrt(s*sqrt(pi));
%     morletWavelet = A * morletWavelet;
%     
%     % Wavelet fft
%     waveletFFT = fft(morletWavelet,Lconv2);
% 
%     % Inverse wavelet fft
%     signalIFT = ifft(waveletFFT .* signalFFT,Lconv2);
%     signalIFT = signalIFT(1:Lconv); % remove the zero-padded signal
%     
%     % Remove the appropriate number of timepoints from the beginning and
%     % end of the time series due to convolution
%     signalIFT = signalIFT(startIndex:endIndex);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%% Normalize Morlet Power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Note: I used this scaling method from Colin Kyle.
    
    % Wavelet fft
    waveletFFT = fft(morletWavelet,Lconv2);

    % Inverse wavelet fft
    signalIFT = ifft((waveletFFT./max(waveletFFT)) .* signalFFT,Lconv2);
    signalIFT = signalIFT(1:Lconv); % remove the zero-padded signal
    
    % Remove the appropriate number of timepoints from the beginning and
    % end of the time series due to convolution
    signalIFT = 2 * signalIFT(startIndex:endIndex);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%     %%%%%%%%%%%%% ANTS Ch 12.6 Scaling Approach %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Note: This doesn't seem to appropriate either, and I don't know why
%     Cohen chose the scaling parameters he did.
%
%     % Wavelet fft
%     waveletFFT = fft(morletWavelet,Lconv2);
% 
%     % Inverse wavelet fft
%     signalIFT = ifft(waveletFFT .* signalFFT,Lconv2);
%     signalIFT = signalIFT(1:Lconv); % remove the zero-padded signal
%     
%     % Remove the appropriate number of timepoints from the beginning and
%     % end of the time series due to convolution and scale with sqrt(s)/10
%     signalIFT = signalIFT(startIndex:endIndex) * sqrt(s)/10;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create power and phase
    pow(i,:) = abs(signalIFT).^2;
    phase(i,:) = atan2(imag(signalIFT),real(signalIFT)); % You can use the function "angle" to calculate this as well
    filtsig(i,:) = real(signalIFT);
    
end

end

