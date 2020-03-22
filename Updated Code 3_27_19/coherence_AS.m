function [coherence,freqs,binSteps] = coherence_AS(LFP1,LFP2,freqRange,timeWindow,stepWindowSize,plotIT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the time-frequency coherence between two signals.
% Thus, it breaks the signal into 2 second windows that are %50 overlapping
% and calculates the coherence during that time. The step windows are
% Hamming windows in order to deal with edge effects from the FFT. Also, it
% handles NaN values. If a NaN is detected within your signal, it analyzes
% the chucks of data surrounding the NaNs. If there are more NaNs than
% signal, then the code places NaNs for the entire signal (in that time
% bin).

%%% Inputs %%%

% LFP1, LFP2: structures with channel data
% Ex: LFP1 = load('/Volumes/SP PHD U3/AmberS/Converted Files/A7_novel1_070116/A7_novel1_070116_07-Mar-2018_CSC12.mat');
% Ex: LFP2 = load('/Volumes/SP PHD U3/AmberS/Converted Files/A7_novel1_070116/A7_novel1_070116_07-Mar-2018_CSC16.mat');
% LFP1 = load('A37_ventral hipp.mat');
% LFP2 = load('Hippocampal EEG.mat');
%
% Ex: freqRange = [2 30];
%
% Size of window of time of interest to detect coherence
% Enter 2 x 2 matrix of the beginning and end of the window of interest for
% each signal. The length of each window must be the same!
% Ex: timeWindow = [0 10; 0 10];
%
% Ex: stepWindowSize = 2;
%
% Ex: plotIT = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input Checks

if LFP1.sFreq ~= LFP2.sFreq
    error('The sampling rates of your signals are different.')
end

timeStamps = LFP1.timestamps - LFP1.timestamps(1);
timeIdx1 = timeStamps >= timeWindow(1,1) & timeStamps < timeWindow(1,2);
timeStamps = LFP2.timestamps - LFP2.timestamps(1);
timeIdx2 = timeStamps >= timeWindow(2,1) & timeStamps < timeWindow(2,2);

if sum(timeIdx1) ~= sum(timeIdx2)
    error('Your singals are not the same length; please adjust your time window matrix.')
end

if timeWindow(1,2) - timeWindow(1,1) < stepWindowSize
    error('Your time window of interest is < 2 seconds making it difficult to get an accurate low frequency estimates.')
end

%% Calculation

x = LFP1.values(timeIdx1);
y = LFP2.values(timeIdx2);
newTimeStamps = 0:1/LFP1.sFreq:(timeWindow(1,2)-timeWindow(1,1))-(1/LFP1.sFreq);

stepPoints = stepWindowSize * LFP1.sFreq;
hamWin = stepPoints/4;
overlap = hamWin/2;
nFFT = [];

binSteps = (0:1000:(length(x) - stepPoints))/1000; % For coherence with smaller steps, replace 1000 (1s) with 250 (quarter of a second).

for iBin = 1:length(binSteps)
    binIdx = newTimeStamps >= binSteps(iBin) & newTimeStamps < binSteps(iBin)+stepWindowSize;
    %[Cxy,F] = mscohereNaN(x,y,stepWindow,overlap,nfft,sFreq)
    [coherence(:,iBin),freqs] = mscohereNaN(x(binIdx),y(binIdx),hamWin,overlap,nFFT,LFP1.sFreq);
end

%% Plot

if plotIT
    
    freqIdx = freqs >= freqRange(1) & freqs <= freqRange(2);
    
    figure
    colormap(jet)
    ax(1) = subplot(211);
    plot(timeStamps(timeIdx1),[x' y'])
    ax(2) = subplot(212);
    imagesc(binSteps,1:sum(freqIdx),coherence(freqIdx,:))
    axis xy
    set(gca,'YTick',1:sum(freqIdx),'YTickLabel',freqs(freqIdx))%,'XTick',1:5:length(binSteps),'XTickLabel',binSteps(1:5:end))
    %xlim([0 120])
    linkaxes(ax,'x')
    
end

%% Extra Notes on Coherence

%%% Alternative way to calculate cxy %%%
% [pxy,~] = cpsd(LFP1.values(timeIdx),LFP2.values(timeIdx),[],[],fqs,1000);
% [pxx,~] = pwelch(LFP1.values(timeIdx),[],[],fqs,LFP1.sFreq);
% [pyy,~] = pwelch(LFP2.values(timeIdx),[],[],fqs,LFP1.sFreq);
% CXY = abs(pxy).^2./(pxx.*pyy);
% isequal(cxy,CXY)

%%% Matt's way to calculate coherece (DAVIS_Coherence_Average.m) %%%
% Use the PDF below to make sense what he is doing
% Notes to calculate coherence: https://atmos.washington.edu/~dennis/552_Notes_6c.pdf
% Briefly:
% Coherence = |Fxy(k)|^2/Fxx*Fyy
% Fxy(k) = Cxk*Cyk*e^(i(thetaxk - thetayk))
% Cxk = sqrt(Fxx(k)) and Cyk = sqrt(Fyy(k))
% power = output from waveletDecomp, tix = indices of time window
% Fxx = nanmean(power(:,tix,1),2); % Power in one
% Fyy = nanmean(power(:,tix,2),2); % Power in two
% specX = abs(nanmean(power(:,tix,1).^0.5 .* power(:,tix,2).^0.5 .* exp(1i*(circ_dist(phase(:,tix,1),phase(:,tix,2)))),2)).^2;
% Cxy = specX./(Fxx.*Fyy);
%
% Matt also applied a median filter, not sure why.
%Cxy = medfilt1(Cxy,5); % The median filter that slides 1 window over and takes the median. Snew = converted power, 20 = window size
%Cxy = medfilt1(Cxy',5)'; % Transposes it

%%% Matlab TF coherence but don't have Wavelet Toolbox %%%
% Using the 30-day trial, it seemed to give basically the same output as above
% [wcoh,wcs,f] = wcoherence(x,y,fs) 
% wcoh = wcoherence(LFP1.values(1:3000),LFP2.values(1:3000),1000);

%%% Unsure about this method %%%
% Time-frequency coherence (https://dsp.stackexchange.com/questions/11503/how-can-i-compute-a-time-frequency-cross-spectrum-in-matlab)
% %fftlen = ;
% overlap = fftlen/2; % 50% overlap
% win = hanning(fftlen);
% X = buffer(x,fftlen,fftlen/2,'nodelay');  % Matrix of overlapping STIs
% numSTIs = size(X,2);
% winX = X.*win(:,ones(1,numSTIs)); % Time-domain windowed STIs
% S = fft(winX,fftlen,1)/fftlen;  % Double-Sided Complex Spectrum Matrix
% %SdB = 20*log10(2*abs(S(1:fftlen/2+1,:)));  % Log Scale Single-Sided Real Spectrum Matrix
% for i = 1:numSTIs
%     Sx1x2(:,:,i) = S1(:,i)'*S2(:,i);
% end

%%% Compare to Chronux Toolbox coherency.m %%%
% x = LFP1.values(1:10000);
% y = LFP2.values(1:10000);
% 
% stepPoints = 2 * LFP1.sFreq;
% hamWin = stepPoints/4;
% overlap = hamWin/2;
% nFFT = [];
% 
% params.tapers = [3 5];
% params.pad = 0;
% params.Fs = LFP1.sFreq;
% params.fpass = [0 50];
% params.err = 0;
% params.trialave = 0;
% 
% binSteps = (0:250:(length(x) - stepPoints))/1000;
% 
% for iBin = 1:length(binSteps)
%     binIdx = newTimeStamps >= binSteps(iBin) & newTimeStamps < binSteps(iBin)+stepWindowSize;
%     %[Cxy,F] = mscohereNaN(x,y,stepWindow,overlap,nfft,sFreq)
%     [Cmatlab(:,iBin),freqs] = mscohereNaN(x(binIdx),y(binIdx),hamWin,overlap,nFFT,LFP1.sFreq);
%     [Cchronux(:,iBin),phi,S12,S1,S2,Fchronux] = coherencyc(x(binIdx),y(binIdx),params);
% end
% 
% figure('Position',[100 100 2000 1000])
% subplot(131)
% imagesc(Cmatlab(1:26,:))
% colorbar; caxis([0 1])
% subplot(132)
% imagesc(Cchronux(1:4:103,:))
% colorbar; caxis([0 1])
% subplot(133)
% imagesc(Cchronux(1:4:103,:)-Cmatlab(1:26,:))
% colorbar; caxis([0 1])

end