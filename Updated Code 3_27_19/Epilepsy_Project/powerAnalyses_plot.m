%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Power Analyses for Gurkoff Lab
%
% Amber Schedlbauer - 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located
%directoryPath = '/Volumes/SP PHD U3/AmberS/Converted Files/NOR/';
%directoryPath = '/Volumes/Gurkoff Lab/Amber/Converted Files/Epileptogenesis/';
directoryPath = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/';

% Enter either:
%   1) Excel file containing the animal information
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
%   2) Mat file containing the "fileInfo" variable (i.e. files2analyze_date.mat)
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/files2analyze_4_2_18.mat';
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% Task to analyze (enter only 1)
session = 'novel1';

% Group to analyze
groupInput = {'SHAM';'PILO'};

% Treatment to analyze
treatmentInput = {'NO'; 'STIM';'BURST'};

% Region to analyze (enter only 1)
electrode = {'LV'};

% Enter a vector of frequencies (log-space) over which you would like to
% calculate power (entire spectrum)
% Example: frequencies = logspace(log10(3),log10(54),24);
frequencies = logspace(log10(2),log10(30),30);

% Enter the frequency band of particular interest
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that you are curious to examine
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [0 120];

% Enter a "1" or "0" if you want to use baseline normalization
% If so, modify the baselineNormalization function to suit your needs.
useBaselineNormalization = 1;

% Path for folder where the results will be saved
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/' session '/' electrode{1} '/'];

% Calculate power for specific behaviors
addSpecificBehaviors = 1;
videoAnalysisFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx';
%behavior = {'NM';'S'};
%behavior = {'M'};
behavior = {'ITO';'IBO'};

%% Identify files to analyze and calculate power

% Grabs the animals of interest
[~,~,ext] = fileparts(codesFile);
if strcmp(ext,'.xlsx')
    %fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
    fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrode,1,[]);
elseif strcmp(ext,'.mat')
    load(codesFile)
end
analysisAnimals = fileInfo(:,1);

% Indices of animals belonging to a treatment group
groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));
groupLabel = cell(1,length(analysisAnimals));
groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; groupLabel(groupIdx{4}) = {'BURST'};

% Read bad channels
channels = readBadChannels(fileInfo,codesFile,2,session);

% Indices of animals belong to a spike group
spikeGroups = {'Low Spike Rate','High Spike Rate'};
spikeThreshold = 2.5;
spikeTable = plotNORspikeTotals(codesFile,'/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/',{'LV'},0,[]);
%     %%% Special Code to add in a spike rate from postpilostim31 for those
%     %%% animals that did not have an NOR recording session
%     keyboard
%     load(['/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/spikeTotals_' analysisAnimals{11} '_postpilostim31_CSC12.mat'],'spikeRate');
%     spikeTable{11,6} = spikeRate;
%     load(['/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/spikeTotals_' analysisAnimals{14} '_postpilostim31_CSC12.mat'],'spikeRate');
%     spikeTable{14,6} = spikeRate;
%     load(['/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/spikeTotals_' analysisAnimals{16} '_postpilostim31_CSC12.mat'],'spikeRate');
%     spikeTable{16,6} = spikeRate;
spikeGroupIdx = cell(1,2);
for iSpike = 1:numel(spikeTable(:,end))
    if spikeTable{iSpike,end} < spikeThreshold % Low spike rate group
        spikeGroupIdx{1} = [spikeGroupIdx{1} iSpike];
    elseif spikeTable{iSpike,end} >= spikeThreshold % High spike rate group
        spikeGroupIdx{2} = [spikeGroupIdx{2} iSpike];
    end
end

% Index of range of frequencies of interest
[~,idxF1] = min(abs(frequencies - frequenciesOfInterest(1)));
[~,idxF2] = min(abs(frequencies - frequenciesOfInterest(2)));

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

% Allocations
powerAll = nan(length(frequencies),length(windowOfInterest(1):windowOfInterest(2)-1)*1000,size(fileInfo,1));
powerAllConvert = nan(length(frequencies),length(windowOfInterest(1):windowOfInterest(2)-1)*1000,size(fileInfo,1));
numWindows = nan(size(fileInfo,1),1);
meanPeakFreqAll = nan(size(fileInfo,1),1);
powerAlign = cell(1,size(fileInfo,1));
powerAlignConvert = nan(length(frequencies),3000,size(fileInfo,1));

%% Loop through animals and removes bad channels and computes power
for iAnimal = 1:size(analysisAnimals)
    
    if ~isempty(channels{iAnimal})
        
        % Load file
        %fileName = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_' session '*' channels{iAnimal} '.mat']);
        fileName = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_*' channels{iAnimal} '.mat']);
        if exist(fileName.name,'file') == 2
            LFP = load(fileName.name); disp(fileName.name);
        else
            continue
        end
        
        % Get time/freq info
        sRate = LFP.sFreq;
        timeStamps = round(LFP.timestamps - LFP.timestamps(1),3); % Grabs total signal time with timepoint(1) = 0
        interestWindowIdx = find(timeStamps >= windowOfInterest(1) & timeStamps < windowOfInterest(2)); % Index of window of interest
        
        % Insert spike information into LFP structure as bad_intervals
        if isfield(LFP,'spike_intervals') % Check to see if this field exists
            LFP.bad_intervals = [LFP.bad_intervals; LFP.spike_intervals];
        end
        
        % Filter raw signal for 60 Hz line noise
        d_60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',sRate);
        LFP.values = filtfilt(d_60,LFP.values);
        %fvtool(d_60) %Look at filter reponse
        
        %%%%%%%%%%%%%%%%%%%%%% POWER CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
        OUT(iAnimal) = power_AS(LFP,frequencies,windowOfInterest,0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Save calculated information for all animals
        powerAll(:,:,iAnimal) = OUT(iAnimal).power(:,interestWindowIdx);

        % Conduct baseline normalization on raw power values
        % [powerDB,powerZ] = baselineNormalization(powerInfo,analysisAnimal,task)
        [powerDB,powerZ] = baselineNormalization(OUT(iAnimal).power,OUT(iAnimal).power(:,interestWindowIdx),analysisAnimals{iAnimal},session);
        powerAllConvert(:,:,iAnimal) = powerDB;
        
        meanPeakFreqAll(iAnimal) = OUT(iAnimal).mean_peak_freq;

        % Sanity Checks: Plot individual raw data and power
        %         figure('Position',[100 100 1500 1000])
        %         ax(1) = subplot(211);
        %         plot(timeStamps(interestWindowIdx),LFP.values(interestWindowIdx))
        %         title('Raw Time Series')
        %         ax(2) = subplot(212);
        %         imagesc(timeStamps(interestWindowIdx),1:length(frequencies),powerAll(:,:,iAnimal))
        %         set(gca,'YTick',1:5:length(frequencies),'YTickLabel',frequencies(1:5:end))
        %         title('Power')
        %         linkaxes(ax,'x')
        %         %saveName = [analysisAnimals{iAnimal} ' ' analysisChan ' Raw Power'];
        %         %saveas(gcf,[saveDirectory saveName])
        
    else
        continue
    end
    
    % Calculates power for the time windows that correspond to the
    % behaviors of interest
    if addSpecificBehaviors
        
        [fileNames,timeWindow] = readNORbehavior(videoAnalysisFile,analysisAnimals(iAnimal),{session},behavior);
        
        % Checks to ensure there are not multiple files
        if isempty(fileNames) || isempty(timeWindow{1})
            continue
        else
            
            timeVec = zeros(1,length(interestWindowIdx));
            
            numWindows(iAnimal) = size(timeWindow{1},1);
            
            for iTime = 1:numWindows(iAnimal)
                
                timeVec(timeWindow{1}(iTime,1)*sRate - 999:timeWindow{1}(iTime,2)*sRate) = 1;
                timeVec = logical(timeVec);
                
                startTime = timeWindow{1}(iTime,1)*sRate - 999;
                endTime = timeWindow{1}(iTime,2)*sRate;
                powerAlign{iAnimal}(:,:,iTime) = powerAll(:,startTime:endTime,iAnimal);
                
            end
            
            % Average trials with-in animal before doing the baseline normalization
            [powerAlignConvert(:,:,iAnimal),~] = baselineNormalization(OUT(iAnimal).power,nanmean(powerAlign{iAnimal},3),analysisAnimals{iAnimal},session);
            
        end
    end
    
    clear fileName LFP fileNames timeWindow tempPower
    
end

%% Logitudinal PSDs
%PSD1 = squeeze(nanmean(powerAll,2));
%powerAll31 = powerAll;
%save('/Users/amberschedlbauer/Desktop/Gurkoff Lab/PowerValues31_raw','powerAll31','groups','groupIdx','frequencies') %'groups','groupIdx'

%% Power Band Ratio
[~,thetaIdx1] = min(abs(frequencies - 6));
[~,thetaIdx2] = min(abs(frequencies - 10));
[~,deltaIdx1] = min(abs(frequencies - 2));
deltaIdx2 = thetaIdx2 - 1;

aveTheta = squeeze(nanmean(nanmean(powerAll(thetaIdx1:thetaIdx2,:,:),2),1));
aveDelta = squeeze(nanmean(nanmean(powerAll(deltaIdx1:deltaIdx2,:,:),2),1));
thetaDeltaRatio = aveTheta./aveDelta;

thetaDelta = [squeeze(nanmean(nanmean(powerAllConvert(thetaIdx1:thetaIdx2,:,:),2),1)) squeeze(nanmean(nanmean(powerAllConvert(deltaIdx1:deltaIdx2,:,:),2),1))];

if addSpecificBehaviors
    aveThetaBehav = nan(1,length(analysisAnimals)); aveDeltaBehav = nan(1,length(analysisAnimals));
    filledIdx = ~cellfun(@isempty,powerAlign); % Need to identify which cells are not empty because cellfun does not play well with empty cells
    aveThetaBehav(filledIdx) = cellfun(@(x) squeeze(nanmean(nanmean(x(thetaIdx1:thetaIdx2,:),2),1)),powerAlign(filledIdx));
    aveDeltaBehav(filledIdx) = cellfun(@(x) squeeze(nanmean(nanmean(x(deltaIdx1:deltaIdx2,:),2),1)),powerAlign(filledIdx));
    thetaDeltaRatioBehav = (aveThetaBehav./aveDeltaBehav)';
    thetaDeltaBehav = [squeeze(nanmean(nanmean(powerAlignConvert(thetaIdx1:thetaIdx2,:,:),2),1)) squeeze(nanmean(nanmean(powerAlignConvert(deltaIdx1:deltaIdx2,:,:),2),1))];

    %aveThetaBehav = squeeze(nanmean(nanmean(powerAlign(thetaIdx1:thetaIdx2,:,:),2),1));
    %aveDeltaBehav = squeeze(nanmean(nanmean(powerAlign(deltaIdx1:deltaIdx2,:,:),2),1));
    %thetaDeltaRatioBehav = aveThetaBehav./aveDeltaBehav;
end

%% Plot

% Set graphing UDIs as empty
UDIs.groups = {};
UDIs.groupColor = {};
UDIs.XAXIS = [];
UDIs.YAXIS = [];
UDIs.COLORMAP = '';
UDIs.XLABEL = {''};
UDIs.YLABEL = '';
UDIs.TITLE = {''};
UDIs.COLORBARTITLE = '';
UDIs.SAVENAME = '';

% Colors for groups
groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};
spikeGroupColor = {[0 0 0];[0.5 0.5 0.5]};

%%% Mean Peak Freq Boxplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 4 Groups %%%
dataCell = {meanPeakFreqAll(groupIdx{1})',meanPeakFreqAll(groupIdx{2})',meanPeakFreqAll(groupIdx{3})',meanPeakFreqAll(groupIdx{4})'};

UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.YLABEL = 'Frequencies (Hz)';

[pval,~,~] = anova1(meanPeakFreqAll,groupLabel,'off'); % Unbalanced 1-way ANOVA for significance
UDIs.TITLE = {['1-Way ANOVA: p = ' num2str(pval)]};

UDIs.SAVENAME = [saveDirectory 'meanPeakFreq_box_allBehaviors_4groups'];

UDIs = plotAnalyses(dataCell,frequencies,UDIs,'boxPlotScatter');
clear data*
%%%

%%% 2 Groups %%%
dataCell = {meanPeakFreqAll(spikeGroupIdx{1})',meanPeakFreqAll(spikeGroupIdx{2})'};

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.YLABEL = 'Frequencies (Hz)';

[~,pval] = ttest2(dataCell{1},dataCell{2}); % Non-paired ttest for significance
UDIs.TITLE = {['2-Sample Ttest: p = ' num2str(pval)]};

UDIs.SAVENAME = [saveDirectory 'meanPeakFreq_box_allBehaviors_2groups'];

UDIs = plotAnalyses(dataCell,frequencies,UDIs,'boxPlotScatter');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Peak Freq Distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[peakPower,peakIdx] = max(powerAllConvert,[],1,'includenan'); % find the maximum of the power in each window
peakIdx(isnan(peakPower)) = NaN;
peakIdx(~isnan(peakIdx));

%%% 4 Groups %%%
dataCell{1} = peakIdx(:,:,groupIdx{1}); dataCell{2} = peakIdx(:,:,groupIdx{2}); dataCell{3} = peakIdx(:,:,groupIdx{3}); dataCell{4} = peakIdx(:,:,groupIdx{4});
UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.XAXIS = frequencies;
UDIs.TITLE = {'Peak Freq Distribution'};
UDIs.SAVENAME = [saveDirectory 'peakFreq_histogram_allBehaviors_4groups'];

UDIs = plotAnalyses(dataCell,frequencies,UDIs,'histogramPlot');
clear data*
%%%

%%% 2 Groups %%%
dataCell{1} = peakIdx(:,:,spikeGroupIdx{1}); dataCell{2} = peakIdx(:,:,spikeGroupIdx{2});
UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.XAXIS = frequencies;
UDIs.TITLE = {'Peak Freq Distribution'};
UDIs.SAVENAME = [saveDirectory 'peakFreq_histogram_allBehaviors_2groups'];

UDIs = plotAnalyses(dataCell,frequencies,UDIs,'histogramPlot');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Theta-Delta Ratio Boxplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 4 Groups %%%
dataCell = {thetaDeltaRatio(groupIdx{1})',thetaDeltaRatio(groupIdx{2})',thetaDeltaRatio(groupIdx{3})',thetaDeltaRatio(groupIdx{4})'};

UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.YLABEL = 'Theta/Delta Ratio';

[pval,~,~] = anova1(thetaDeltaRatio,groupLabel,'off'); % Unbalanced 1-way ANOVA for significance
UDIs.TITLE = {['1-Way ANOVA: p = ' num2str(pval)]};

UDIs.SAVENAME = [saveDirectory 'thetaDeltaRatio_box_allBehaviors_4groups'];

UDIs = plotAnalyses(dataCell,frequencies,UDIs,'boxPlotScatter');
clear data*
%%%

%%% 2 Groups %%%
dataCell = {thetaDeltaRatio(spikeGroupIdx{1})',thetaDeltaRatio(spikeGroupIdx{2})'};

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.YLABEL = 'Theta/Delta Ratio';

[~,pval] = ttest2(dataCell{1},dataCell{2}); % Non-paired ttest for significance
UDIs.TITLE = {['2-Sample Ttest: p = ' num2str(pval)]};

UDIs.SAVENAME = [saveDirectory 'thetaDeltaRatio_box_allBehaviors_2groups'];

UDIs = plotAnalyses(dataCell,frequencies,UDIs,'boxPlotScatter');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Average spectrogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 4 Groups %%%
data = cat(3,nanmean(powerAllConvert(:,:,groupIdx{1}),3),nanmean(powerAllConvert(:,:,groupIdx{2}),3),nanmean(powerAllConvert(:,:,groupIdx{3}),3),nanmean(powerAllConvert(:,:,groupIdx{4}),3));
dataCell{1} = data(:,:,2) - data(:,:,1); dataCell{2} = data(:,:,3) - data(:,:,1); dataCell{3} = data(:,:,4) - data(:,:,1);

UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.XAXIS = timeStamps(interestWindowIdx);
UDIs.YAXIS = frequencies;
UDIs.COLORMAP = 'jet';
UDIs.XLABEL = {'','Time (s)',''};
UDIs.YLABEL = 'Frequencies (Hz)';
UDIs.TITLE = {'PILO-SHAM';'STIM-SHAM';'BURST-SHAM'};
UDIs.COLORBARTITLE = 'dB';
UDIs.SAVENAME = [saveDirectory 'spectrogram_allBehaviors_4groups'];

UDIs = plotAnalyses(dataCell,frequencies,UDIs,'matrix');
clear data*
%%%

%%% 2 Groups %%%%
data = cat(3,nanmean(powerAllConvert(:,:,spikeGroupIdx{1}),3),nanmean(powerAllConvert(:,:,spikeGroupIdx{2}),3));
dataCell{1} = data(:,:,2) - data(:,:,1);

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.XAXIS = timeStamps(interestWindowIdx);
UDIs.YAXIS = frequencies;
UDIs.COLORMAP = 'jet';
UDIs.XLABEL = {'Time (s)'};
UDIs.YLABEL = 'Frequencies (Hz)';
UDIs.TITLE = {[spikeGroups{2} ' - ' spikeGroups{1}]};
UDIs.COLORBARTITLE = 'dB';
UDIs.SAVENAME = [saveDirectory 'spectrogram_allBehaviors_2groups'];

UDIs = plotAnalyses(dataCell,frequencies,UDIs,'matrix');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Mean PSD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSD = squeeze(nanmean(powerAllConvert,2));

%%% 4 groups %%%
dataCell{1} = nanmean(PSD(:,groupIdx{1}),2); dataCell{2} = nanmean(PSD(:,groupIdx{2}),2); dataCell{3} = nanmean(PSD(:,groupIdx{3}),2); dataCell{4} = nanmean(PSD(:,groupIdx{4}),2);

UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.XAXIS = frequencies;
UDIs.XLABEL = {'Frequencies (Hz)'};
UDIs.YLABEL = 'dB';
UDIs.SAVENAME = [saveDirectory 'PSDmean_allBehaviors_4groups'];

UDIs = plotAnalyses(dataCell,frequencies,UDIs,'PSD');
clear data*
%%%

%%% 2 groups %%%
dataCell{1} = PSD(:,spikeGroupIdx{1})'; dataCell{2} = PSD(:,spikeGroupIdx{2})';

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.XAXIS = frequencies;
UDIs.XLABEL = {'Frequencies (Hz)'};
UDIs.YLABEL = 'dB';
UDIs.SAVENAME = [saveDirectory 'PSDmean_allBehaviors_2groups'];

UDIs = plotAnalyses(dataCell,frequencies,UDIs,'PSDsig');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Separate Theta/Delta Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[100 100 500 500])
hold on
h1 = plot(1:2,thetaDelta(groupIdx{1},:),'-','Color',[groupColor{1} 0.5],'LineWidth',2);
h2 = plot(1:2,nanmean(thetaDelta(groupIdx{1},:)),'-','Color','k','LineWidth',4);
h3 = plot(3:4,thetaDelta(groupIdx{2},:)','-','Color',[groupColor{2} 0.5],'LineWidth',2);
h4 = plot(3:4,nanmean(thetaDelta(groupIdx{2},:)),'-','Color','k','LineWidth',4);
h5 = plot(5:6,thetaDelta(groupIdx{1},:),'-','Color',[groupColor{3} 0.5],'LineWidth',2);
h6 = plot(5:6,nanmean(thetaDelta(groupIdx{1},:)),'-','Color','k','LineWidth',4);
h7 = plot(7:8,thetaDelta(groupIdx{2},:)','-','Color',[groupColor{4} 0.5],'LineWidth',2);
h8 = plot(7:8,nanmean(thetaDelta(groupIdx{2},:)),'-','Color','k','LineWidth',4);
xlim([0 9])
ylim([-7 2])
ylabel('Power (dB)')
set(gca,'XTick',1:8,'XTickLabel',{'\Theta','\Delta','\Theta','\Delta','\Theta','\Delta','\Theta','\Delta'})
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
legend([h1(1) h3(1) h5(1) h7(1)],groups)
save2pdf([saveDirectory 'thetaDelta_allBehaviors_4groups'])

figure('Position',[100 100 500 500])
hold on
h1 = plot(1:2,thetaDelta(spikeGroupIdx{1},:),'-','Color',[spikeGroupColor{1} 0.5],'LineWidth',2);
h2 = plot(1:2,nanmean(thetaDelta(spikeGroupIdx{1},:)),'-','Color','r','LineWidth',4);
h3 = plot(3:4,thetaDelta(spikeGroupIdx{2},:)','-','Color',[spikeGroupColor{2} 0.5],'LineWidth',2);
h4 = plot(3:4,nanmean(thetaDelta(spikeGroupIdx{2},:)),'-','Color','r','LineWidth',4);
xlim([0 5])
ylim([-7 2])
ylabel('Power (dB)')
set(gca,'XTick',1:4,'XTickLabel',{'\Theta','\Delta','\Theta','\Delta'})
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
legend([h1(1) h3(1)],spikeGroups)
save2pdf([saveDirectory 'thetaDelta_allBehaviors_2groups'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if addSpecificBehaviors
    
    %%% Theta-Delta Ratio Boxplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 4 Groups %%%
    dataCell = {thetaDeltaRatioBehav(groupIdx{1})',thetaDeltaRatioBehav(groupIdx{2})',thetaDeltaRatioBehav(groupIdx{3})',thetaDeltaRatioBehav(groupIdx{4})'};
    
    UDIs.groups = groups;
    UDIs.groupColor = groupColor;
    UDIs.YLABEL = 'Theta/Delta Ratio';
    
    [pval,~,~] = anova1(thetaDeltaRatioBehav,groupLabel,'off'); % Unbalanced 1-way ANOVA for significance
    UDIs.TITLE = {['1-Way ANOVA: p = ' num2str(pval)]};
    
    UDIs.SAVENAME = [saveDirectory 'thetaDeltaRatio_box_4groups'];
    
    UDIs = plotAnalyses(dataCell,frequencies,UDIs,'boxPlotScatter');
    clear data*
    %%%
    
    %%% 2 Groups %%%
    dataCell = {thetaDeltaRatioBehav(spikeGroupIdx{1})',thetaDeltaRatioBehav(spikeGroupIdx{2})'};
    
    UDIs.groups = spikeGroups;
    UDIs.groupColor = spikeGroupColor;
    UDIs.YLABEL = 'Theta/Delta Ratio';
    
    [~,pval] = ttest2(dataCell{1},dataCell{2}); % Non-paired ttest for significance
    UDIs.TITLE = {['2-Sample Ttest: p = ' num2str(pval)]};
    
    UDIs.SAVENAME = [saveDirectory 'thetaDeltaRatio_box_2groups'];
    
    UDIs = plotAnalyses(dataCell,frequencies,UDIs,'boxPlotScatter');
    clear data*
    %%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Average spectrogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 4 Groups %%%
    data = cat(3,nanmean(powerAlignConvert(:,:,groupIdx{1}),3),nanmean(powerAlignConvert(:,:,groupIdx{2}),3),nanmean(powerAlignConvert(:,:,groupIdx{3}),3),nanmean(powerAlignConvert(:,:,groupIdx{4}),3));
    dataCell{1} = data(:,:,2) - data(:,:,1); dataCell{2} = data(:,:,3) - data(:,:,1); dataCell{3} = data(:,:,4) - data(:,:,1);
    
    UDIs.groups = groups;
    UDIs.groupColor = groupColor;
    UDIs.XAXIS = -0.999:0.001:2;
    UDIs.YAXIS = frequencies;
    UDIs.COLORMAP = 'jet';
    UDIs.XLABEL = {'','Time (s)',''};
    UDIs.YLABEL = 'Frequencies (Hz)';
    UDIs.TITLE = {'PILO-SHAM';'STIM-SHAM';'BURST-SHAM'};
    UDIs.COLORBARTITLE = 'dB';
    UDIs.SAVENAME = [saveDirectory 'spectrogram_4groups'];
    
    UDIs = plotAnalyses(dataCell,frequencies,UDIs,'matrix');
    clear data*
    %%%
    
    %%% 2 Groups %%%%
    data = cat(3,nanmean(powerAlignConvert(:,:,spikeGroupIdx{1}),3),nanmean(powerAlignConvert(:,:,spikeGroupIdx{2}),3));
    dataCell{1} = data(:,:,2) - data(:,:,1);
    
    UDIs.groups = spikeGroups;
    UDIs.groupColor = spikeGroupColor;
    UDIs.XAXIS = -0.999:0.001:2;
    UDIs.YAXIS = frequencies;
    UDIs.COLORMAP = 'jet';
    UDIs.XLABEL = {'Time (s)'};
    UDIs.YLABEL = 'Frequencies (Hz)';
    UDIs.TITLE = {[spikeGroups{2} ' - ' spikeGroups{1}]};
    UDIs.COLORBARTITLE = 'dB';
    UDIs.SAVENAME = [saveDirectory 'spectrogram_2groups'];
    
    UDIs = plotAnalyses(dataCell,frequencies,UDIs,'matrix');
    clear data*
    %%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Mean PSD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PSD = squeeze(nanmean(powerAlignConvert,2));
    
    %%% 4 groups %%%
    dataCell{1} = nanmean(PSD(:,groupIdx{1}),2); dataCell{2} = nanmean(PSD(:,groupIdx{2}),2); dataCell{3} = nanmean(PSD(:,groupIdx{3}),2); dataCell{4} = nanmean(PSD(:,groupIdx{4}),2);
    
    UDIs.groups = groups;
    UDIs.groupColor = groupColor;
    UDIs.XAXIS = frequencies;
    UDIs.XLABEL = {'Frequencies (Hz)'};
    UDIs.YLABEL = 'dB';
    UDIs.SAVENAME = [saveDirectory 'PSDmean_4groups'];
    
    UDIs = plotAnalyses(dataCell,frequencies,UDIs,'PSD');
    clear data*
    %%%
    
    %%% 2 groups %%%
    dataCell{1} = PSD(:,spikeGroupIdx{1})'; dataCell{2} = PSD(:,spikeGroupIdx{2})';
    
    UDIs.groups = spikeGroups;
    UDIs.groupColor = spikeGroupColor;
    UDIs.XAXIS = frequencies;
    UDIs.XLABEL = {'Frequencies (Hz)'};
    UDIs.YLABEL = 'dB';
    UDIs.SAVENAME = [saveDirectory 'PSDmean_2groups'];
    
    UDIs = plotAnalyses(dataCell,frequencies,UDIs,'PSDsig');
    clear data*
    %%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Separate Theta/Delta Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure('Position',[100 100 500 500])
    hold on
    h1 = plot(1:2,thetaDeltaBehav(groupIdx{1},:),'-','Color',[groupColor{1} 0.5],'LineWidth',2);
    h2 = plot(1:2,nanmean(thetaDeltaBehav(groupIdx{1},:)),'-','Color','k','LineWidth',4);
    h3 = plot(3:4,thetaDeltaBehav(groupIdx{2},:)','-','Color',[groupColor{2} 0.5],'LineWidth',2);
    h4 = plot(3:4,nanmean(thetaDeltaBehav(groupIdx{2},:)),'-','Color','k','LineWidth',4);
    h5 = plot(5:6,thetaDeltaBehav(groupIdx{1},:),'-','Color',[groupColor{3} 0.5],'LineWidth',2);
    h6 = plot(5:6,nanmean(thetaDeltaBehav(groupIdx{1},:)),'-','Color','k','LineWidth',4);
    h7 = plot(7:8,thetaDeltaBehav(groupIdx{2},:)','-','Color',[groupColor{4} 0.5],'LineWidth',2);
    h8 = plot(7:8,nanmean(thetaDeltaBehav(groupIdx{2},:)),'-','Color','k','LineWidth',4);
    xlim([0 9])
    %ylim([-7 2])
    ylabel('Power (dB)')
    set(gca,'XTick',1:8,'XTickLabel',{'\Theta','\Delta','\Theta','\Delta','\Theta','\Delta','\Theta','\Delta'})
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    legend([h1(1) h3(1) h5(1) h7(1)],groups)
    save2pdf([saveDirectory 'thetaDelta_4groups'])
    
    figure('Position',[100 100 500 500])
    hold on
    h1 = plot(1:2,thetaDeltaBehav(spikeGroupIdx{1},:),'-','Color',[spikeGroupColor{1} 0.5],'LineWidth',2);
    h2 = plot(1:2,nanmean(thetaDeltaBehav(spikeGroupIdx{1},:)),'-','Color','r','LineWidth',4);
    h3 = plot(3:4,thetaDeltaBehav(spikeGroupIdx{2},:)','-','Color',[spikeGroupColor{2} 0.5],'LineWidth',2);
    h4 = plot(3:4,nanmean(thetaDeltaBehav(spikeGroupIdx{2},:)),'-','Color','r','LineWidth',4);
    xlim([0 5])
    %ylim([-7 2])
    ylabel('Power (dB)')
    set(gca,'XTick',1:4,'XTickLabel',{'\Theta','\Delta','\Theta','\Delta'})
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    legend([h1(1) h3(1)],spikeGroups)
    save2pdf([saveDirectory 'thetaDelta_2groups'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end

soundAlarm

%% MISC EXTRA

%%% STATS %%%
%dataCell = {thetaDeltaRatioBehav(groupIdx{1})',thetaDeltaRatioBehav(groupIdx{2})',thetaDeltaRatioBehav(groupIdx{4})'}; %,thetaDeltaRatioBehav(groupIdx{3})'
%dataCell = {thetaDeltaRatioBehav(spikeGroupIdx{1})',thetaDeltaRatioBehav(spikeGroupIdx{2})'};

% maxSize = max(cellfun(@numel,dataCell));    % Get the maximum vector size
% fcn = @(x) [x nan(1,maxSize-numel(x))];  % Create an anonymous function
% rmat = cellfun(fcn,dataCell,'UniformOutput',false);  % Pad each cell with NaNs
% rmat = vertcat(rmat{:})';
%
% [p,tbl,stats] = anova1(rmat);

%%% PLOTS %%%
% 
% %%% Mean + Individual PSD for 4 groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD = PSD';
% figure('Position',[100 100 1500 1000])
% for iGroup = 1:4
%     subplot(2,2,iGroup)
%     plot(frequencies,PSD(:,groupIdx{iGroup}),'Color',groupColor{iGroup})
%     hold on
%     plot(frequencies,nanmean(PSD(:,groupIdx{iGroup}),2),'LineWidth',6,'Color',groupColor{iGroup})
%     set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
%     xlabel('Frequencies (Hz)')
%     ylabel('dB')
%     title(groups{iGroup})
%     ylim([-10 5])
%     vline(frequenciesOfInterest,'k--')
% end
% save2pdf([saveDirectory 'PSDindividual_4groups'])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Mean + Individual PSD for 2 groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('Position',[100 100 1500 500])
% for iGroup = 1:2
%     subplot(1,2,iGroup)
%     plot(frequencies,PSD(:,spikeGroupIdx{iGroup}),'Color',spikeGroupColor{iGroup})
%     hold on
%     plot(frequencies,nanmean(PSD(:,spikeGroupIdx{iGroup}),2),'LineWidth',6,'Color',spikeGroupColor{iGroup})
%     set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
%     xlabel('Frequencies (Hz)')
%     ylabel('dB')
%     title(spikeGroups{iGroup})
%     ylim([-10 5])
%     vline(frequenciesOfInterest,'k--')
% end
% save2pdf([saveDirectory 'PSDindividual_2groups'])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
