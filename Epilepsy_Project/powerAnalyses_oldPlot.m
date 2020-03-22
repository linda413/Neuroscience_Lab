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
directoryPath = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/';
%directoryPath = '/Volumes/SP PHD U3/AmberS/Converted Files/NOR/';

% Enter either:
%   1) Excel file containing the animal information
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
%   2) Mat file containing the "fileInfo" variable (i.e. files2analyze_date.mat)
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/files2analyze_4_2_18.mat';
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% Task to analyze (enter only 1)
session = 'novel2';

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

% Indices of animals belong to a spike group
spikeGroups = {'Low Spike Rate';'High Spike Rate'};
spikeThreshold = 2.5;
spikeTable = plotNORspikeTotals(codesFile,'/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/',{'LV'},0,[]);
spikeGroupIdx = cell(1,2);
for iSpike = 1:numel(spikeTable(:,end))
    if spikeTable{iSpike,end} < spikeThreshold % Low spike rate group
        spikeGroupIdx{1} = [spikeGroupIdx{1} iSpike];
    elseif spikeTable{iSpike,end} >= spikeThreshold % High spike rate group
        spikeGroupIdx{2} = [spikeGroupIdx{2} iSpike];
    end
end

% Read bad channels
channels = readBadChannels(fileInfo,codesFile,2,session);

% Index of range of frequencies of interest
[~,idxF1] = min(abs(frequencies - frequenciesOfInterest(1)));
[~,idxF2] = min(abs(frequencies - frequenciesOfInterest(2)));

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
        fileName = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_' session '*' channels{iAnimal} '.mat']);
        if exist(fileName.name,'file') == 2
            LFP = load(fileName.name); disp(fileName.name);
        else
            continue
        end
        
        % Get time/freq info
        sRate = LFP.sFreq;
        timeStamps = LFP.timestamps - LFP.timestamps(1); % Grabs total signal time with timepoint(1) = 0
        interestWindowIdx = find(timeStamps > windowOfInterest(1) & timeStamps <= windowOfInterest(2)); % Index of window of interest
        
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
        [powerDB,powerZ(:,:,iAnimal)] = baselineNormalization(OUT(iAnimal).power,OUT(iAnimal).power(:,interestWindowIdx),analysisAnimals{iAnimal},session);
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

%% Power Band Ratio
[~,thetaIdx1] = min(abs(frequencies - 6));
[~,thetaIdx2] = min(abs(frequencies - 10));
[~,deltaIdx1] = min(abs(frequencies - 2));
deltaIdx2 = thetaIdx2 - 1;

aveTheta = squeeze(nanmean(nanmean(powerAll(thetaIdx1:thetaIdx2,:,:),2),1));
aveDelta = squeeze(nanmean(nanmean(powerAll(deltaIdx1:deltaIdx2,:,:),2),1));
thetaDeltaRatio = aveTheta./aveDelta;

if addSpecificBehaviors
    aveThetaBehav = nan(1,length(analysisAnimals)); aveDeltaBehav = nan(1,length(analysisAnimals));
    filledIdx = ~cellfun(@isempty,powerAlign); % Need to identify which cells are not empty because cellfun does not play well with empty cells
    aveThetaBehav(filledIdx) = cellfun(@(x) squeeze(nanmean(nanmean(x(thetaIdx1:thetaIdx2,:),2),1)),powerAlign(filledIdx));
    aveDeltaBehav(filledIdx) = cellfun(@(x) squeeze(nanmean(nanmean(x(deltaIdx1:deltaIdx2,:),2),1)),powerAlign(filledIdx));
    thetaDeltaRatioBehav = (aveThetaBehav./aveDeltaBehav)';
    %aveThetaBehav = squeeze(nanmean(nanmean(powerAlign(thetaIdx1:thetaIdx2,:,:),2),1));
    %aveDeltaBehav = squeeze(nanmean(nanmean(powerAlign(deltaIdx1:deltaIdx2,:,:),2),1));
    %thetaDeltaRatioBehav = aveThetaBehav./aveDeltaBehav;
end

%% Plot

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

% Colors for treatment groups
groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};
spikeGroupColor = {[0 0 0];[0.5 0.5 0.5]};

UDIs.groupColor = groupColor;
UDIs.XAXIS = [];
UDIs.YAXIS = [];
UDIs.COLORMAP = [];
UDIs.XLABEL = '';
UDIs.YLABEL = '';
UDIs.TITLE = '';
UDIs.COLORBARTITLE = '';
UDIs.SAVENAME = '';

%%% Mean Peak Freq Boxplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 4 Groups %%%
dataCell = {meanPeakFreqAll(groupIdx{1})',meanPeakFreqAll(groupIdx{2})',meanPeakFreqAll(groupIdx{3})',meanPeakFreqAll(groupIdx{4})'};

figure('Position',[100 100 500 500])
plotBoxScatter(dataCell,groups,groupColor)
set(gca,'Position',[0.13 0.12 0.8 0.8])
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
ylabel('Frequencies (Hz)')

% Unbalanced 1-way ANOVA for significance
groupLabel = cell(1,length(analysisAnimals));
groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; groupLabel(groupIdx{4}) = {'BURST'};
[pval,~,~] = anova1(meanPeakFreqAll,groupLabel,'off');
title(['1-Way ANOVA: p = ' num2str(pval)])

save2pdf([saveDirectory 'meanPeakFreq_box_4groups'])
%%%

%%% 2 Groups %%%
figure('Position',[100 100 500 500])
dataCell = {meanPeakFreqAll(spikeGroupIdx{1})',meanPeakFreqAll(spikeGroupIdx{2})'};
plotBoxScatter(dataCell,spikeGroups,spikeGroupColor)
set(gca,'Position',[0.13 0.12 0.8 0.8])
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
ylabel('Frequencies (Hz)')

% Non-paired ttest for significance
[~,pval] = ttest2(dataCell{1},dataCell{2});
%text(1,2,['2-Sample Ttest: p = ' num2str(pval)])
%title(['Mean Peak Freq - ' session ' - ' electrode{1}])
title(['2-Sample Ttest: p = ' num2str(pval)])

save2pdf([saveDirectory 'meanPeakFreq_box_2groups'])
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Peak Freq Distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[peakPower,peakIdx] = max(powerAllConvert,[],1,'includenan'); % find the maximum of the power in each window
peakIdx(isnan(peakPower)) = NaN;
figure('Position',[100 100 2000 500])
hold on
for iGroup = 1:4
    temp = peakIdx(:,:,groupIdx{iGroup});
    %subplot(4,1,iGroup)
    histogram(temp(:),'FaceColor',groupColor{iGroup},'FaceAlpha',0.5) 
end
set(gca,'XTick',1:length(frequencies),'XTickLabel',round(frequencies,1))
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
legend(groups)
save2pdf([saveDirectory 'peakFreq_histogram_4groups'])

figure('Position',[100 100 2000 500])
hold on
for iGroup = 1:2
    temp = peakIdx(:,:,spikeGroupIdx{iGroup});
    %subplot(4,1,iGroup)
    histogram(temp(:),'FaceColor',spikeGroupColor{iGroup},'FaceAlpha',0.5) 
end
set(gca,'XTick',1:length(frequencies),'XTickLabel',round(frequencies,1))
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
legend(spikeGroups)
save2pdf([saveDirectory 'peakFreq_histogram_2groups'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Theta-Delta Ratio Boxplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Position',[100 100 1500 500])

%%% 4 Groups %%%
figure('Position',[100 100 500 500])
%subplot(1,2,1)
dataCell = {thetaDeltaRatio(groupIdx{1})',thetaDeltaRatio(groupIdx{2})',thetaDeltaRatio(groupIdx{3})',thetaDeltaRatio(groupIdx{4})'};
plotBoxScatter(dataCell,groups,groupColor)
set(gca,'Position',[0.13 0.12 0.8 0.8])
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
ylabel('Theta/Delta Ratio')

% Unbalanced 1-way ANOVA for significance
groupLabel = cell(1,length(analysisAnimals));
groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; groupLabel(groupIdx{4}) = {'BURST'};
[pval,~,~] = anova1(thetaDeltaRatio,groupLabel,'off');
%text(1,2,['1-Way ANOVA: p = ' num2str(pval)])
title(['1-Way ANOVA: p = ' num2str(pval)])

save2pdf([saveDirectory 'thetaDeltaRatio_box_allBehaviors_4groups'])
%%%

%%% 2 Groups %%%
figure('Position',[100 100 500 500])
%subplot(1,2,2)
dataCell = {thetaDeltaRatio(spikeGroupIdx{1})',thetaDeltaRatio(spikeGroupIdx{2})'};
plotBoxScatter(dataCell,spikeGroups,spikeGroupColor)
set(gca,'Position',[0.13 0.12 0.8 0.8])
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
ylabel('Theta/Delta Ratio')

% Non-paired ttest for significance
[~,pval] = ttest2(dataCell{1},dataCell{2});
%text(1,2,['2-Sample Ttest: p = ' num2str(pval)])
title(['2-Sample Ttest: p = ' num2str(pval)])

save2pdf([saveDirectory 'thetaDeltaRatio_box_allBehaviors_2groups'])
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Average spectrogram for 4 groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = cat(3,nanmean(powerAllConvert(:,:,groupIdx{1}),3),nanmean(powerAllConvert(:,:,groupIdx{2}),3),nanmean(powerAllConvert(:,:,groupIdx{3}),3),nanmean(powerAllConvert(:,:,groupIdx{4}),3));
%dataCell{1} = data(:,:,2) - data(:,:,1); dataCell{2} = data(:,:,3) - data(:,:,1); dataCell{3} = data(:,:,4) - data(:,:,1);
dataDiff(:,:,1) = data(:,:,2) - data(:,:,1); dataDiff(:,:,2) = data(:,:,3) - data(:,:,1); dataDiff(:,:,3) = data(:,:,4) - data(:,:,1);
colorbarLims = [min(dataDiff(:)) max(dataDiff(:))];

figure('Position',[100 100 2000 500])
for iGroup = 1:size(dataDiff,3)
    ax(iGroup) = subplot(1,size(dataDiff,3),iGroup);
    imagesc(timeStamps(interestWindowIdx),1:length(frequencies),dataDiff(:,:,iGroup))
    axis xy
    set(gca,'YTick',1:5:length(frequencies),'YTickLabel',[])
    caxis(colorbarLims)
    if iGroup == 1
        ax(iGroup).Position = [0.13 0.12 0.2 0.8];
        set(gca,'YTickLabel',round(frequencies(1:5:end),1))
        ylabel('Frequencies (Hz)')
    elseif iGroup == 2
        ax(iGroup).Position = [0.34 0.12 0.2 0.8];
        xlabel('Time (s)')
    elseif iGroup == 3
        cH = colorbar;
        title(cH,'dB')
        ax(iGroup).Position = [0.55 0.12 0.2 0.8];
    end
    hline([idxF1 idxF2],{'--k','--k'})
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    title([groups{iGroup+1} ' - ' groups{1}])
end
save2pdf([saveDirectory 'spectrogram_allBehaviors_4groups'])
clear data dataDiff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Average spectrogram for 2 groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = cat(3,nanmean(powerAllConvert(:,:,spikeGroupIdx{1}),3),nanmean(powerAllConvert(:,:,spikeGroupIdx{2}),3));
% colorbarLims = [min(data(:)) max(data(:))];
% 
% figure('Position',[100 100 1500 500])
% for iGroup = 1:length(spikeGroupIdx)
%     ax2(iGroup) = subplot(1,2,iGroup);
%     imagesc(timeStamps(interestWindowIdx),1:length(frequencies),data(:,:,iGroup))
%     axis xy
%     set(gca,'YTick',1:5:length(frequencies),'YTickLabel',frequencies(1:5:end))
%     set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
%     hline([idxF1 idxF2],{'--k','--k'})
%     cH = colorbar;
%     cH.Label.String = 'dB';
%     caxis(colorbarLims)
%     xlabel('Time (s)')
%     ylabel('Frequencies (Hz)')
%     title(spikeGroups{iGroup})
% end
% linkaxes(ax2)

dataDiff = data(:,:,2) - data(:,:,1);
figure('Position',[100 100 565 500])
imagesc(timeStamps(interestWindowIdx),1:length(frequencies),dataDiff)
axis xy
hline([idxF1 idxF2],{'--k','--k'})
set(gca,'YTick',1:5:length(frequencies),'YTickLabel',round(frequencies(1:5:end),1))
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
cH = colorbar;
title(cH,'dB')
caxis(colorbarLims)
set(gca,'Position',[0.13 0.12 0.708 0.8])
xlabel('Time (s)')
ylabel('Frequencies (Hz)')
title([spikeGroups{2} ' - ' spikeGroups{1}])
save2pdf([saveDirectory 'spectrogram_allBehaviors_2groups'])
clear data dataDiff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Mean PSD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure('Position',[100 100 1500 500])

%%% Mean PSD for 4 groups %%%
figure('Position',[100 100 500 500])
%subplot(121)
PSD = squeeze(nanmean(powerAllConvert,2));
for iGroup = 1:4
    plot(frequencies,nanmean(PSD(:,groupIdx{iGroup}),2),'LineWidth',6,'Color',groupColor{iGroup})
    hold on
end
vline(frequenciesOfInterest,'k--')
legend(groups)
set(gca,'Position',[0.13 0.12 0.8 0.8])
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
xlabel('Frequencies (Hz)')
ylabel('dB')
save2pdf([saveDirectory 'PSDmean_allBehaviors_4groups'])
%%%

%%% Mean PSD for 2 groups %%%

% USE FOR NO SIGNIFICANCE
% % figure('Position',[100 100 750 500])
% % %subplot(122)
% % PSD = squeeze(nanmean(powerAllConvert,2));
% % for iGroup = 1:2
% %     plot(frequencies,nanmean(PSD(:,spikeGroupIdx{iGroup}),2),'LineWidth',6,'Color',spikeGroupColor{iGroup})
% %     hold on
% % end
% % vline(frequenciesOfInterest,'k--')
% % legend(spikeGroups)
% % set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
% % xlabel('Frequencies (Hz)')
% % ylabel('dB')
% % save2pdf([saveDirectory 'PSDmean_allBehaviors_2groups'])

% USE FOR SIGNIFICANCE
% Uses a Welch's ttest for only the frequencies in the theta range and
% corrects for multiple comparisons using Bonferroni
PSD = squeeze(nanmean(powerAllConvert,2))';
[h,p] = ttest2(PSD(spikeGroupIdx{1},:),PSD(spikeGroupIdx{2},:),'Vartype','unequal','Alpha',0.05/numel(idxF1:idxF2));
freqWin = zeros(1,length(frequencies)); freqWin(idxF1:idxF2) = 1; h = logical(h.*freqWin);

dataMeans(:,1) = nanmean(PSD(spikeGroupIdx{1},:));
dataMeans(:,2) = nanmean(PSD(spikeGroupIdx{2},:));

figure('Position',[100 100 500 500]);
fill([frequencies(h) flip(frequencies(h))],[dataMeans(h,2)' flip(dataMeans(h,1))'],'k','LineStyle','none')
alpha(0.15)
hold on
hL = plot(frequencies,dataMeans,'LineWidth',4);
hL(1).Color = spikeGroupColor{1}; hL(2).Color = spikeGroupColor{2};
vline(frequenciesOfInterest,'k--')
legend(hL,spikeGroups)
set(gca,'Position',[0.13 0.12 0.8 0.8])
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
xlabel('Frequencies (Hz)')
ylabel('dB')
save2pdf([saveDirectory 'PSDmean_allBehaviors_2groups'])
%%%

%save2pdf([saveDirectory 'PSDmean_' session '_' electrode{1} '_baselineNormal'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if addSpecificBehaviors
    
    behaviorList = cell(1,length(behavior));
    behaviorList(1:end-1) = {' '};
    behaviorList = [behavior behaviorList']';
    
    %%% Theta-Delta Ratio Boxplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %figure('Position',[100 100 1500 500])
    
    %%% 4 Groups %%%
    figure('Position',[100 100 500 500])
    %subplot(1,2,1)
    dataCell = {thetaDeltaRatioBehav(groupIdx{1})',thetaDeltaRatioBehav(groupIdx{2})',thetaDeltaRatioBehav(groupIdx{3})',thetaDeltaRatioBehav(groupIdx{4})'};
    plotBoxScatter(dataCell,groups,groupColor)
    set(gca,'Position',[0.13 0.12 0.8 0.8])
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    %title(['Theta/Delta Ratio Behavior: ' [behaviorList{:}] ' - ' session ' - ' electrode{1}])
    ylabel('Theta/Delta Ratio')
    
    % Unbalanced 1-way ANOVA for significance
    groupLabel = cell(1,length(analysisAnimals));
    groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; groupLabel(groupIdx{4}) = {'BURST'};
    [pval,~,~] = anova1(thetaDeltaRatioBehav,groupLabel,'off');
    %text(1,2,['1-Way ANOVA: p = ' num2str(pval)])
    title(['1-Way ANOVA: p = ' num2str(pval)])
    %[results,means] = multcompare(stats,'CType','bonferroni')
    
    save2pdf([saveDirectory 'thetaDeltaRatio_box_4groups'])
    %%%
    
    %%% 2 Groups %%%
    figure('Position',[100 100 500 500])
    %subplot(1,2,2)
    dataCell = {thetaDeltaRatioBehav(spikeGroupIdx{1})',thetaDeltaRatioBehav(spikeGroupIdx{2})'};
    plotBoxScatter(dataCell,spikeGroups,spikeGroupColor)
    set(gca,'Position',[0.13 0.12 0.8 0.8])
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    %title(['Theta/Delta Ratio Behavior: ' [behaviorList{:}] ' - ' session ' - ' electrode{1}])
    ylabel('Theta/Delta Ratio')
    
    % Non-paired ttest for significance
    [h,pval] = ttest2(dataCell{1},dataCell{2});
    %text(1,2,['2-Sample Ttest: p = ' num2str(pval)])
    title(['2-Sample Ttest: p = ' num2str(pval)])
    
    save2pdf([saveDirectory 'thetaDeltaRatio_box_2groups'])
    %%%
    
    %save2pdf([saveDirectory 'thetaDeltaRatio_box'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Average spectrogram for 4 groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure('Position',[100 100 3000 300])
%     %data = cat(3,nanmean(powerConvertAlign(:,:,groupIdx{1}),3),nanmean(powerConvertAlign(:,:,groupIdx{2}),3),nanmean(powerConvertAlign(:,:,groupIdx{3}),3),nanmean(powerConvertAlign(:,:,groupIdx{4}),3));
%     colorbarLims = [min(data(:)) max(data(:))];
%     
%     for iGroup = 1:length(groupIdx)
%         ax1(iGroup) = subplot(1,4,iGroup);
%         imagesc(-0.999:0.001:2,1:length(frequencies),nanmean(powerAlignConvert(:,:,groupIdx{iGroup}),3))
%         axis xy
%         set(gca,'YTick',1:5:length(frequencies),'YTickLabel',round(frequencies(1:5:end),1))
%         set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
%         hline([idxF1 idxF2],{'--k','--k'})
%         cH = colorbar;
%         cH.Label.String = 'dB';
%         caxis(colorbarLims)
%         xlabel('Time (s)')
%         ylabel('Frequency (Hz)')
%         title(groups{iGroup})
%         %saveas(gcf,[saveDirectory 'spectrogram_' groups{iGroup}],'jpeg')
%     end
%     linkaxes(ax1)
%     save2pdf([saveDirectory 'spectrogram_4groups'])
    
    data = cat(3,nanmean(powerAlignConvert(:,:,groupIdx{1}),3),nanmean(powerAlignConvert(:,:,groupIdx{2}),3),nanmean(powerAlignConvert(:,:,groupIdx{3}),3),nanmean(powerAlignConvert(:,:,groupIdx{4}),3));
    dataDiff(:,:,1) = data(:,:,2) - data(:,:,1); dataDiff(:,:,2) = data(:,:,3) - data(:,:,1); dataDiff(:,:,3) = data(:,:,4) - data(:,:,1); 
    colorbarLims = [min(dataDiff(:)) max(dataDiff(:))];

    figure('Position',[100 100 2000 500])
    for iGroup = 1:size(dataDiff,3)
        ax(iGroup) = subplot(1,size(dataDiff,3),iGroup);
        imagesc(-0.999:0.001:2,1:length(frequencies),dataDiff(:,:,iGroup))
        axis xy
        set(gca,'YTick',1:5:length(frequencies),'YTickLabel',[])
        caxis(colorbarLims)
        if iGroup == 1
            ax(iGroup).Position = [0.13 0.12 0.2 0.8];
            set(gca,'YTickLabel',round(frequencies(1:5:end),1))
            ylabel('Frequencies (Hz)')
        elseif iGroup == 2
            ax(iGroup).Position = [0.34 0.12 0.2 0.8];
            xlabel('Time (s)')
        elseif iGroup == 3
            cH = colorbar;
            title(cH,'dB')
            ax(iGroup).Position = [0.55 0.12 0.2 0.8];
        end
        hline([idxF1 idxF2],{'--k','--k'})
        set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
        title([groups{iGroup+1} ' - ' groups{1}])
    end
    save2pdf([saveDirectory 'spectrogram_4groups'])
    clear data dataDiff
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Average spectrogram for 2 groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure('Position',[100 100 1500 500])
%     data = cat(3,nanmean(powerAlignConvert(:,:,spikeGroupIdx{1}),3),nanmean(powerAlignConvert(:,:,spikeGroupIdx{2}),3));
%     colorbarLims = [min(data(:)) max(data(:))];
%     
%     for iGroup = 1:length(spikeGroupIdx)
%         ax2(iGroup) = subplot(1,2,iGroup);
%         imagesc(-1:0.001:1.999,1:length(frequencies),data(:,:,iGroup))
%         axis xy
%         set(gca,'YTick',1:5:length(frequencies),'YTickLabel',round(frequencies(1:5:end),1))
%         set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
%         hline([idxF1 idxF2],{'--k','--k'})
%         cH = colorbar;
%         cH.Label.String = 'dB';
%         caxis(colorbarLims)
%         xlabel('Time (s)')
%         ylabel('Frequencies (Hz)')
%         title(spikeGroups{iGroup})
%         %saveas(gcf,[saveDirectory 'spectrogram_' spikeGroups{iGroup}],'jpeg')
%     end
%     linkaxes(ax2)
%     save2pdf([saveDirectory 'spectrogram_2groups'])
    
    data = cat(3,nanmean(powerAlignConvert(:,:,spikeGroupIdx{1}),3),nanmean(powerAlignConvert(:,:,spikeGroupIdx{2}),3));
    dataDiff = data(:,:,2) - data(:,:,1);
    figure('Position',[100 100 565 500])
    imagesc(-1:0.001:1.999,1:length(frequencies),dataDiff)
    axis xy
    hline([idxF1 idxF2],{'--k','--k'})
    set(gca,'YTick',1:5:length(frequencies),'YTickLabel',round(frequencies(1:5:end),1))
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    cH = colorbar;
    title(cH,'dB')
    caxis(colorbarLims)
    set(gca,'Position',[0.13 0.12 0.708 0.8])
    xlabel('Time (s)')
    ylabel('Frequencies (Hz)')
    title([spikeGroups{2} ' - ' spikeGroups{1}])
    save2pdf([saveDirectory 'spectrogram_2groups'])
    clear data dataDiff
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Mean PSD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %figure('Position',[100 100 1500 500])
    
    %%% Mean PSD for 4 groups %%%
    figure('Position',[100 100 500 500])
    %subplot(121)
    PSD = squeeze(nanmean(powerAlignConvert,2));
    %PSD = squeeze(nanmean(powerAll,2));
    for iGroup = 1:4
        plot(frequencies,nanmean(PSD(:,groupIdx{iGroup}),2),'LineWidth',6,'Color',groupColor{iGroup})
        hold on
    end
    vline(frequenciesOfInterest,'k--')
    legend(groups)
    set(gca,'Position',[0.13 0.12 0.8 0.8])
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    xlabel('Frequencies (Hz)')
    ylabel('dB')
    save2pdf([saveDirectory 'PSDmean_4groups'])
    %%%
    
    %%% Mean PSD for 2 groups %%%
    
    % USE FOR SIGNIFICANCE
    % Uses a Welch's ttest for only the frequencies in the theta range and
    % corrects for multiple comparisons using Bonferroni
    PSD = squeeze(nanmean(powerAlignConvert,2))';
    [h,p] = ttest2(PSD(spikeGroupIdx{1},:),PSD(spikeGroupIdx{2},:),'Vartype','unequal','Alpha',0.05/numel(idxF1:idxF2));
    freqWin = zeros(1,length(frequencies)); freqWin(idxF1:idxF2) = 1; h = logical(h.*freqWin);
    
    dataMeans(:,1) = nanmean(PSD(spikeGroupIdx{1},:));
    dataMeans(:,2) = nanmean(PSD(spikeGroupIdx{2},:));
    
    figure('Position',[100 100 500 500]);
    if sum(h) == 1
        plot([frequencies(h) frequencies(h)],[dataMeans(h,2) flip(dataMeans(h,1))],'Color',[0.5 0.5 0.5 0.15],'LineWidth',5);
    else
        fill([frequencies(h) flip(frequencies(h))],[dataMeans(h,2)' flip(dataMeans(h,1))'],'k','LineStyle','none')
    end
    alpha(0.15)
    hold on
    hL = plot(frequencies,dataMeans,'LineWidth',4);
    hL(1).Color = spikeGroupColor{1}; hL(2).Color = spikeGroupColor{2};
    vline(frequenciesOfInterest,'k--')
    legend(hL,spikeGroups)
    set(gca,'Position',[0.13 0.12 0.8 0.8])
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    xlabel('Frequencies (Hz)')
    ylabel('dB')
    save2pdf([saveDirectory 'PSDmean_2groups'])
    %%%
    
    %save2pdf([saveDirectory 'PSDmean_' session '_' electrode{1} '_baselineNormal'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Mean + Individual PSD for 4 groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PSD = PSD';
    figure('Position',[100 100 1500 1000])
    for iGroup = 1:4
        subplot(2,2,iGroup)
        plot(frequencies,PSD(:,groupIdx{iGroup}),'Color',groupColor{iGroup})
        hold on
        plot(frequencies,nanmean(PSD(:,groupIdx{iGroup}),2),'LineWidth',6,'Color',groupColor{iGroup})
        set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
        xlabel('Frequencies (Hz)')
        ylabel('dB')
        title(groups{iGroup})
        ylim([-10 5])
        vline(frequenciesOfInterest,'k--')
    end
    save2pdf([saveDirectory 'PSDindividual_4groups'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Mean + Individual PSD for 2 groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Position',[100 100 1500 500])
    for iGroup = 1:2
        subplot(1,2,iGroup)
        plot(frequencies,PSD(:,spikeGroupIdx{iGroup}),'Color',spikeGroupColor{iGroup})
        hold on
        plot(frequencies,nanmean(PSD(:,spikeGroupIdx{iGroup}),2),'LineWidth',6,'Color',spikeGroupColor{iGroup})
        set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
        xlabel('Frequencies (Hz)')
        ylabel('dB')
        title(spikeGroups{iGroup})
        ylim([-10 5])
        vline(frequenciesOfInterest,'k--')
    end
    save2pdf([saveDirectory 'PSDindividual_2groups'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

soundAlarm

%% Stats

%dataCell = {thetaDeltaRatioBehav(groupIdx{1})',thetaDeltaRatioBehav(groupIdx{2})',thetaDeltaRatioBehav(groupIdx{4})'}; %,thetaDeltaRatioBehav(groupIdx{3})'
%dataCell = {thetaDeltaRatioBehav(spikeGroupIdx{1})',thetaDeltaRatioBehav(spikeGroupIdx{2})'};

% maxSize = max(cellfun(@numel,dataCell));    % Get the maximum vector size
% fcn = @(x) [x nan(1,maxSize-numel(x))];  % Create an anonymous function
% rmat = cellfun(fcn,dataCell,'UniformOutput',false);  % Pad each cell with NaNs
% rmat = vertcat(rmat{:})';
%
% [p,tbl,stats] = anova1(rmat);
