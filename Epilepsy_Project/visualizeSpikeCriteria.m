% Visualize spike criteria
% Run countInterictalSpikes_AS.m and keyboard in Animal Loop

hFig = figure('Position',[100 100 1500 500]);

% Raw data and standard deviation threshold
%ts = LFP.timestamps-LFP.timestamps(1);
ax(1) = subplot(411);
plot(values,'k') % plot 1st to set axes correctly
ylim([min(values(:))-0.25 max(values(:))+0.25])
hold on

if ~isempty(posLocs)
    vline(posLocs,'--g')
end

if ~isempty(negLocs)
    vline(negLocs,'--g')
end

% Power criteria
normPower = powerNaN(:,:)';
normPower = bsxfun(@rdivide,normPower,max(normPower));
aveNormPower = mean(normPower,2);
spikeThresh = quantile(aveNormPower,numQuantile);

[~,powLocs,~,~] = findpeaks(double(aveNormPower),sFreq,'MinPeakHeight',spikeThresh,'MinPeakDistance',1);
powLocs1 = round(powLocs*sFreq);

ax(2) = subplot(412);
plot(values,'k') % plot 1st to set axes correctly
hold on
if ~isempty(powLocs1)
    vline(powLocs1,'--b')
end

% Closest power peaks to std peaks
combLocs = sort([posLocs negLocs]);
if ~isempty(combLocs)
    for iClosest = 1:length(powLocs1)
        [minDist,closestIdx] = min(abs(combLocs - powLocs1(iClosest)));
        if minDist <= 50
            powLocs2(iClosest) = combLocs(closestIdx);
        end
    end
end

ax(3) = subplot(413);
plot(values,'k') % plot 1st to set axes correctly
hold on
if ~isempty(powLocs2)
    vline(powLocs2,'--c')
end

% Slow wave complex criteria
dLowPass = designfilt('lowpassfir','PassbandFrequency',5,'StopbandFrequency',10,'PassbandRipple',1,'StopbandAttenuation',60,'SampleRate',1000);
lowPassData = filtfilt(dLowPass,values);
lowPassDataNaN = lowPassData;
lowPassDataNaN(artifactReject) = NaN;
spikeThresh = numSTDslow * nanstd(lowPassDataNaN);
[posPeaksLow,posLocsLow,~,~] = findpeaks(lowPassData,sFreq,'MinPeakHeight',spikeThresh,'MinPeakDistance',1);
[negPeaksLow,negLocsLow,~,~] = findpeaks(lowPassData*-1,sFreq,'MinPeakHeight',spikeThresh,'MinPeakDistance',1);

posLocsLow = round(posLocsLow*sFreq); negLocsLow = round(negLocsLow*sFreq);

ax(4) = subplot(414);
plot(values,'k') % plot 1st to set axes correctly
hold on

if ~isempty(posLocsLow)
    vline(posLocs,'--r')
end

if ~isempty(negLocsLow)
    vline(negLocs,'--r')
end

linkaxes(ax)