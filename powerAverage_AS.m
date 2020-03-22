function OUT = powerAverage_AS(LFP,freqs,timeWindow)

% OUT = powerAverage_AS(LFP,freqs,timeWindow)
% LFP: LFP structure or a file
% freqs: vector of frequencies of interest
% timeWindow: time window of interest in seconds
%   Ex: [startTime endTime]
%   NOTE: If you have n windows of interest, please provide a n x 2 matrix
%   of start and end times of windows.

%% Check inputs
if isempty(LFP)
    LFP = uigetfile('Select LFP');
    LFP = load(LFP);
end

if nargin < 3
    timeWindow = [0 (LFP.timestamps(end) - LFP.timestamps(1))];
end

%% Extract power of signal in LFP

waveletNum = 6;

% Handle bad data by interpolation (better handles edge effects from fft)
[interpolatedData,badIdx] = interpolateBadData(LFP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NO EDGE BUFFER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract power from raw LFP data using wavelet decomposition over
% range of frequencies specifiefd by freqs
%[~,power,filteredSignal] = waveletDecomp_AS(freqs,interpolatedData,LFP.sFreq,waveletNum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NO EDGE BUFFER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDGE BUFFER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Buffer the EEG data to ensure no edges artifacts at beginning and end
[bufferedData,bufferBins] = addMirroredBuffers(interpolatedData,freqs,LFP.sFreq);

% Extract power from raw LFP data using wavelet decomposition over
% range of frequencies specifiefd by freqs
[~,power,~] = waveletDecomp_AS(freqs,bufferedData,LFP.sFreq,waveletNum);    

% Trim the buffers
power = power(:,(bufferBins + 1):(length(power) - bufferBins));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDGE BUFFER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NaN artifact time points
power(:,badIdx) = NaN;
badVec = isnan(power(1,:));

%% Average power over time windows

OUT = zeros(size(timeWindow,1),1);

for iTime = 1:size(timeWindow,1)
    
    ts = LFP.timestamps - LFP.timestamps(1); % make timestamps relative to the start of recording
    tix = ts >= timeWindow(iTime,1) & ts < timeWindow(iTime,2); % binary vector of all timepoints in this window
    
    % Warning if more than half of the signal is artifact
    if mean(badVec(tix)) > 0.5
        display(['WARNING: Average in window ' num2str(iTime) ' is more that 50% rejected!']);
    end
    
    OUT(iTime,1) = nanmean(nanmean(power(:,tix),2)); %Average over time and powerband
    
    display(['Average power in ' num2str(freqs(1)) 'Hz to ' num2str(freqs(end)) 'Hz is ' num2str(OUT(iTime,1)) ' mV^2'])
    
end
    
end

