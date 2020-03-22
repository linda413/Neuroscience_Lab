function [Pepisode,powerDurWindow,power] = pEpisode_AS(LFP,freqs,interestFreqs,backWindow,interestWindow,saveDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%   LFP: structure that contains the data to analyze, output from
%   DAVIS_Manual_Artifact_Reject.m
%   freqs: Vector of frequencies over which you would like to perform
%   pEpisode
%   interestFreqs: [minFreq maxFreq] -> band of frequencies of particular
%   interest
%   backWindow: [minTime maxTime] -> window of time to calculate background
%   spectrum
%   interestWindow: [minTime maxTime] -> window of time to calculate
%   pEpisode
%   saveDir: directory of where to save the output and plots
% OUTPUTS
%   Pepisode: vector of percet time of window of interest spent in the
%   different frequencies outlined in background frequencies
%   powerDurWindow: binary matrix of time points and frequencies that meet
%   the power and duration thresholds
%   power: power output from the spectral decomposition
%   Figures: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters to adjust if needed
waveletNum = 6; % number of cycles for the wavelet for spectral decomposition
percentile = 0.95; %percentile to threshold the chi-square (2) distribution 
durationCycles = 3; % number of consecutive cycles that must exceed the power threshold

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
[~,power,filteredSignal] = waveletDecomp_AS(freqs,bufferedData,LFP.sFreq,waveletNum);    

% Trim the buffers
power = power(:,(bufferBins + 1):(length(power) - bufferBins));
filteredSignal = filteredSignal(:,(bufferBins + 1):(length(filteredSignal) - bufferBins));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDGE BUFFER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NaN artifact time points
power(:,badIdx) = NaN;

% Find the average power for each frequency over the background window
timeStamps = LFP.timestamps - LFP.timestamps(1); %grabs total signal time with timepoint(1) = 0
timeStamps = round(timeStamps,3); % AS added 12_7_18
backWindowIdx = timeStamps >= backWindow(1) & timeStamps < backWindow(2);
powerMean = nanmean(power(:,backWindowIdx),2); % average over the specified window

%%%%%%%%%%%%%%%%%%%%%%%%%% Chi-square calculation%%%%%%%%%%%%%%%%%%%%%%%%%%
powerLog = log10(powerMean); % transform to log space
pv = polyfit(log10(freqs),powerLog',1); % linear regression fit -> pv is [slope intercept]
fitMeans = 10.^(polyval(pv,log10(freqs))); % fitted mean powers in non-log space

% Rationale (Matthew Schmit) behind calculation on line below:
% The distribution of powers is shaped like a chisquare(2) distribution, 
% but it's actual mean and variance are different at different frequencies.
% So, we need a conversion factor:
%   p = X(2)*c
%   where p = power distribution
%         X(2) = chisquare distribution
%         c = conversion factor
% The mean of the converted distribution is the fitted power, and the mean
% of the X(2) distribution is 2. To convert bewteen the two means, the 
% equation looks like this:
%   p_fitted = mean(X(2))*c
%   p_fitted = 2*c
%   p_fitted/2 = c
% Thus, the converstion factor from X(2) to powers:
%   p = X(2)*p_fitted/2
% This equation would let us convert from X(2) values to power values.
% The chi2inv then gets the X(2) value at the percentile of interest and 
% multiplying by the fitted power divided by two converts from X(2) to
% power. 

% The code calculates chi2inv(0:0.001:0.999,2) which finds the value for 
% all probabilites in an inverse chi-square(2) cdf. For a chi-square
% distribution with 2 degrees of freedom, the pdf = 1/2*e^(-x/2). 
% The cdf = 1 - e^(-x/2). The inverse of the cdf function is
% cdf_inv = -2*ln(1-x). For example, chi2inv(0.95,2) and -2*log(1-0.95)
% give the same value.

nBins = 1000;
thresholdAll = chi2inv(0:(1/nBins):(1-(1/nBins)),2)' * (fitMeans/2);
%%%%%%%%%%%%%%%%%%%%%%%%%% Chi-square calculation%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finds the powers that exceed the power threshold
thresholdSingle = repmat(thresholdAll(nBins*percentile,:)',1,size(power,2));
powerThresh = power > thresholdSingle; % binary vector indicating which timepoints and frequencies met the calculated power threshold

% Finds the segments that meet the duration threshold
for i = 1:length(freqs)
   
   foo = bwlabel(powerThresh(i,:)); % finds the connected components in the powerThresh matrix at a particular frequency
   
   if sum(foo) ~= 0
       counts = histcounts(foo,'BinMethod','integers'); % finds the size of the connected components
       counts = counts(2:end); % remove counts of zero
       idx = find(counts/LFP.sFreq < durationCycles/freqs(i)); % finds which components do not exceed the duration threshold
       
       for j = idx
           foo(foo == j) = 0; % sets the too small components to zero
       end
   end
   
   durationThresh(i,:) = foo;
   
   clear foo counts idx
      
end

durationThresh = durationThresh > 0;

% Finds the above threshold sements for the window of interest
interestWindowIdx = find(timeStamps > interestWindow(1) & timeStamps <= interestWindow(2)); %finds the timepoints within the window of interest
powerDurWindow = durationThresh(:,interestWindowIdx);

% Finds p-episode in the time window of interest
Pepisode = nanmean(powerDurWindow,2)*100;

%%%%%%%%%%%%%%%%%%%%%%%%% Ploting and Saving Output %%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(saveDir)
    
    if ~exist(saveDir,'dir')
        mkdir(saveDir)
    end
    
    %%% Plots power and duration segements %%%
    
    % segments that meet the power threshold
    figure('Position',[100 100 1500 500])
    colormap(hot)
    subplot(131)
    imagesc(powerThresh(:,backWindowIdx))
    set(gca,'XTick',1:sum(backWindowIdx)-1:sum(backWindowIdx),'XTickLabel',backWindow,'YTick',1:length(freqs),'YTickLabel',freqs)
    xlabel('Time (sec)')
    ylabel('Frequencies (Hz)')
    title('Signal that Reaches Power Threshold (Entire Backgroud Window)')
    
    % segments that meets both power and duration thresholds
    subplot(132)
    imagesc(durationThresh(:,backWindowIdx))
    set(gca,'XTick',1:sum(backWindowIdx)-1:sum(backWindowIdx),'XTickLabel',backWindow,'YTick',1:length(freqs),'YTickLabel',freqs)
    xlabel('Time (sec)')
    ylabel('Frequencies (Hz)')
    title('Signal that Reaches Power and Duration Threshold (Entire Backgroud Window)')
    
    % segments that meets both power and duration thresholds over window of interest
    subplot(133)
    imagesc(powerDurWindow)
    set(gca,'XTick',1:size(powerDurWindow,2)-1:size(powerDurWindow,2),'XTickLabel',interestWindow,'YTick',1:length(freqs),'YTickLabel',freqs)
    xlabel('Time (sec)')
    ylabel('Frequencies (Hz)')
    title('Signal that Reaches Power and Duration Threshold (Window of Interest)')
    
    save2pdf([saveDir LFP.Description '_PowerDurationThresholds'])
    close(gcf)
    
    %%% Plots p-episode %%%
    figure
    plot(freqs,Pepisode)
    xlabel('Frequencies (Hz)')
    ylabel('Percent Time')
    title('P-Episode (Window of Interest)')
    save2pdf([saveDir LFP.Description '_Pepisode_' num2str(interestWindow(1)) '_' num2str(interestWindow(2)) 'sec'])
    close(gcf)
    
    %%% Plot output in range of frequencies of interest %%%
    [~,idxF1] = min(abs(freqs - interestFreqs(1)));
    [~,idxF2] = min(abs(freqs - interestFreqs(2)));
    
    ctr = 1;
    
    for iPlot = idxF1:idxF2
        
        threshVal = sqrt(thresholdSingle(iPlot,1));
        
        figure(ctr+2)
        shadeH = fill([0 length(interestWindowIdx) length(interestWindowIdx) 0],[threshVal threshVal -1*threshVal -1*threshVal],'r');
        shadeH.FaceAlpha = 0.25;
        hold on
        plot(filteredSignal(iPlot,interestWindowIdx))
        set(gca,'XTick',1:size(powerDurWindow,2)-1:size(powerDurWindow,2),'XTickLabel',interestWindow)
        xlabel('Time (s)')
        %set(gca,'XTick',0:5000:length(interestWindowIdx),'XTickLabel',(0:5000:length(interestWindowIdx))/1000)
        saveas(gcf,[saveDir LFP.Description '_FilteredSignal_' num2str(freqs(iPlot)) 'Hz.fig'])
        close(gcf)
        
        ctr = ctr + 1;
    end
    
end

end