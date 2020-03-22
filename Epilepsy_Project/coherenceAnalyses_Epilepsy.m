%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coherence Analyses for Gurkoff Lab

% Amber Schedlbauer - 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located
directoryPath = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/';
%directoryPath = '/Volumes/SP PHD U3/AmberS/Converted Files/NOR/';

% Animal information excel file
codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';

% Task to analyze (enter only 1)
session = 'novel2';

% Group to analyze
groupInput = {'SHAM';'PILO'};

% Treatment to analyze
treatmentInput = {'NO'; 'STIM';'BURST'};

% Regions to analyze (enter only 2)
electrodes = {'LV','LPFC'};

% Enter the frequency range for coherence anlaysis
freqRange = [2 30];

% Enter the frequency band of particular interest
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that you are curious to examine
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [0 120];

% Path for folder where the results will be saved
%saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/'];
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/' session '/' electrodes{1} '/'];

% Calculate power for specific behaviors
addSpecificBehaviors = 1;
videoAnalysisFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx';
%behavior = {'NM';'S'};
%behavior = {'M'};
behavior = {'ITO';'IBO'};

%% Identify files to analyze and calculate coherence

% Grabs the animals of interest
%fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrodes{1},1,[]);
channels1 = readBadChannels(fileInfo,codesFile,2,session);
fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrodes{2},1,[]);
channels2 = readBadChannels(fileInfo,codesFile,2,session);
analysisAnimals = fileInfo(:,1);

% Indices of animals belonging to a treatment group
groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));

% Indices of animals belong to a spike group (Include only for epilepsy)
spikeGroups = {'Low Spike Rate';'High Spike Rate'};
spikeTable = plotNORspikeTotals(codesFile,'/Users/amberschedlbauer/Desktop/Gurkoff Lab/Spikes/',{'LV'},0,[]);

% Spike grouping based on cluster
spikeClusterGroup = [2 2 2 1 2 1 1 NaN 1 2 NaN 2 1 NaN NaN NaN 2 1 1 1 2 1 2 2 2 1 1 NaN 2 1 2 1 2 1 1 2 1 1 NaN 2 2 NaN 2 2 1 NaN 1 2 2 1];
spikeGroupIdx{1} = find(spikeClusterGroup == 1);
spikeGroupIdx{2} = find(spikeClusterGroup == 2);

% Spike grouping based on rate
%spikeThreshold = 2.5;
% for iSpike = 1:numel(spikeTable(:,end))
%     if spikeTable{iSpike,end} < spikeThreshold % Low spike rate group
%         spikeGroupIdx{1} = [spikeGroupIdx{1} iSpike];
%     elseif spikeTable{iSpike,end} >= spikeThreshold % High spike rate group
%         spikeGroupIdx{2} = [spikeGroupIdx{2} iSpike];
%     end
% end

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

%% Loop through animals to calculate coherence

ctr = 1;

for iAnimal = 1:length(analysisAnimals)
    
    if ~isempty(channels1{iAnimal}) && ~isempty(channels2{iAnimal}) % Skips animals if there is a bad channel
        
        % Load FILES
        fileName1 = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_' session '*' channels1{iAnimal} '.mat']);
        LFP1 = load(fileName1.name);
        
        fileName2 = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_' session '*' channels2{iAnimal} '.mat']);
        LFP2 = load(fileName2.name);
        
        % Get time/freq info
        sRate = LFP1.sFreq;
        timeStamps = round(LFP1.timestamps - LFP1.timestamps(1),3);
        interestWindowIdx = find(timeStamps == windowOfInterest(1)):(find(timeStamps == windowOfInterest(2))-1);
        
        % Insert spike information into LFP structure as bad_intervals
        if isfield(LFP1,'spike_intervals') % Check to see if this field exists
            LFP1.bad_intervals = [LFP1.bad_intervals; LFP1.spike_intervals];
        end
        
        if isfield(LFP2,'spike_intervals') % Check to see if this field exists
            LFP2.bad_intervals = [LFP2.bad_intervals; LFP2.spike_intervals];
        end
        
        % Filter raw signal for 60 Hz line noise
        d_60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',sRate);
        LFP1.values = filtfilt(d_60,LFP1.values);
        LFP2.values = filtfilt(d_60,LFP2.values);
        %fvtool(d_60) %Look at filter reponse
        
        %%%%%%%%%%%%%%%%%%%% COHERENCE CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%
        [interpolatedData1,badIdx1] = interpolateBadData(LFP1); LFP1.values = interpolatedData1;
        [interpolatedData2,badIdx2] = interpolateBadData(LFP2); LFP2.values = interpolatedData2;
        
        %[coherence,freqs] = coherence_AS(LFP1,LFP2,freqRange,timeWindow,stepWindowSize,plotIT)
        [coherence,freqs,timeBins] = coherence_AS(LFP1,LFP2,freqRange,[windowOfInterest; windowOfInterest],4,0);
        timeBins = timeBins + 2; % Step window size = 4 (see above), so the timing of the coherence estimate is shifted
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Frequencies of interest
        freqIdx = freqs >= freqRange(1) & freqs <= freqRange(2);
        freqThetaIdx = freqs >= frequenciesOfInterest(1) & freqs <= frequenciesOfInterest(2);
        
        % Find timeBins with artifact and remove
        artifactReject = zeros(1,length(interpolatedData1));
        artifactReject(unique([badIdx1 badIdx2])) = 1;
        artifactReject = logical(artifactReject(interestWindowIdx));
        newTime = 0:1/LFP1.sFreq:(windowOfInterest(2) - windowOfInterest(1))-(1/LFP1.sFreq);
        artifactTimes = round(newTime(artifactReject),3);
        for iReshape = 1:length(timeBins)
            timeBinsReshaped(iReshape,:) = timeBins(iReshape):(1/LFP1.sFreq):(timeBins(iReshape) + 1 - 1/LFP1.sFreq);
        end
        timeBinsReshaped = round(timeBinsReshaped,3);
        removeIdx = find(sum(ismember(timeBinsReshaped,artifactTimes),2));
        coherence(:,removeIdx) = NaN;
        
        % Save calculated information for all animals
        if ctr == 1
            coherenceAll = nan([size(coherence) length(analysisAnimals)]);
            behaviorMatAll = nan([sum(freqIdx) size(coherence,2) length(analysisAnimals)]);
            coherenceAlignAll = nan([sum(freqIdx) 9 length(analysisAnimals)]);
            ctr = 2;
        end
        coherenceAll(:,:,iAnimal) = coherence;
        
        if addSpecificBehaviors
            
            % Read in video behavior and add to x labels for "Sanity Check" below %%%
            [~,txt,~] = xlsread(videoAnalysisFile,1,'C3:DS102');
            videoIdx = find(cellfun(@(x) ~isempty(x),strfind(txt(:,1),analysisAnimals{iAnimal})) .* cellfun(@(x) ~isempty(x),strfind(txt(:,1),session)));
            txt = txt(:,2:end);
            for iLabel = 1:length(timeBins)
                if  mod(timeBins(iLabel),1) == 0
                    newLabels{iLabel} = [num2str(timeBins(iLabel)) ' ' txt{videoIdx(1),round(timeBins(iLabel))}];
                else
                    newLabels{iLabel} = num2str(timeBins(iLabel));
                end
            end
            
            %%% One way to visualize the coherence for each animal
            behaviorMat = nan(sum(freqIdx),length(timeBins));
            for iMat = 1:length(timeBins)
                if strcmp(txt{videoIdx(1),timeBins(iMat)},'NO DATA')
                    behaviorMat(:,iMat) = NaN;
                elseif strcmp(txt{videoIdx(1),timeBins(iMat)},'M')
                    behaviorMat(:,iMat) = 1;
                elseif strcmp(txt{videoIdx(1),timeBins(iMat)},'NM')
                    behaviorMat(:,iMat) = 2;
                elseif strcmp(txt{videoIdx(1),timeBins(iMat)},'S')
                    behaviorMat(:,iMat) = 3;
                elseif strcmp(txt{videoIdx(1),timeBins(iMat)},'SNOT') || strcmp(txt{videoIdx(1),timeBins(iMat)},'SNOB')
                    behaviorMat(:,iMat) = 4;
                elseif strcmp(txt{videoIdx(1),timeBins(iMat)},'ITO') || strcmp(txt{videoIdx(1),timeBins(iMat)},'IBO')
                    behaviorMat(:,iMat) = 5;
                elseif strcmp(txt{videoIdx(1),timeBins(iMat)},'ATO') || strcmp(txt{videoIdx(1),timeBins(iMat)},'ABO')
                    behaviorMat(:,iMat) = 6;
                elseif strcmp(txt{videoIdx(1),timeBins(iMat)},'LTO') || strcmp(txt{videoIdx(1),timeBins(iMat)},'LBO')
                    behaviorMat(:,iMat) = 7;
                end
            end
            
            behaviorMatAll(:,:,iAnimal) = behaviorMat;
            
            numWin = find(behaviorMat(1,:) == 1);
            behaviorIdx = [ 1 (find(diff(numWin) > 4)+1)];
            for i = 1:length(behaviorIdx)
                tempIdx = (behaviorIdx(i)-4):(behaviorIdx(i)+4);
                if sum(tempIdx < 1) || sum(tempIdx > size(behaviorMat,2))
                    coherenceAlign(:,:,i) = nan(sum(freqIdx),numel(tempIdx));
                else
                    coherenceAlign(:,:,i) = coherence(freqIdx,tempIdx);
                end
            end
            coherenceAlignAll(:,:,iAnimal) = nanmean(coherenceAlign,3);
            clear coherenceAlign
            
            % View counts of certain behaviors 
            figure
            x = coherence(9,:); %x(x<0.4) = NaN;
            y = behaviorMat(6,:);
            nanIdx = isnan(x + y);
            x(nanIdx) = []; y(nanIdx) = [];
            hHan = histogram2(x(:),y(:),0.05:0.1:1,0.5:1:8,'DisplayStyle','tile','ShowEmptyBins','on');
            set(gca,'YTick',1:7,'YTickLabel',{'M';'NM';'S';'SNOT/SNOB';'ITO/IBO';'ATO/ABO';'LTO/LBO'})
            colorbar
            title(analysisAnimals{iAnimal})
            histCounts(:,:,iAnimal) = hHan.Values';
            close(gcf)
            
        end
        
    else
        continue
    end
    
    % Sanity Checks: Plot individual raw data and coherence
    %         figure('Position',[100 100 1500 1000])
    %         colormap(jet)
    %         ax(1) = subplot(211);
    %         plot(timeStamps(interestWindowIdx),[interpolatedData1(interestWindowIdx)' interpolatedData2(interestWindowIdx)'])
    %         set(gca,'XTick',timeBins,'XTickLabel',timeBins)
    %         legend(electrodes)
    %         title(analysisAnimals{iAnimal})
    %         ax(2) = subplot(212);
    %         imagesc(timeBins+0.50,1:sum(freqIdx),coherence(freqIdx,:)) % Added the 0.250 so that the coherence for the first bin would be plotted at the tick mark and not have the tick mark in the middle of the bin
    %         axis xy
    %         set(gca,'YTick',1:sum(freqIdx),'YTickLabel',freqs(freqIdx),'XTick',timeBins,'XTickLabel',newLabels)
    %         %set(gca,'YTick',1:sum(freqIdx),'YTickLabel',freqs(freqIdx),'XTick',timeBins,'XTickLabel',timeBins)
    %         linkaxes(ax,'x')
    %         waitfor(gcf)
    
    clear fileName* LFP* interpolatedData* coherence artifactReject timeBinsReshaped
    
end

%% Plot

% Colors for treatment groups
groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};
spikeGroupColor = {[0 0 0];[0.5 0.5 0.5]};

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

%%% Matrix of coherence before and after specific behaviors %%%%%%%%%%%%%%%

%%% 2 Groups %%%
dataCell{1} = nanmean(coherenceAlignAll(:,:,spikeGroupIdx{1}),3);
dataCell{2} = nanmean(coherenceAlignAll(:,:,spikeGroupIdx{2}),3);

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.XAXIS = -4:4;
UDIs.YAXIS = freqs(freqIdx);
UDIs.COLORMAP = 'jet';
UDIs.XLABEL = {'Time (s)';'Time (s)'};
UDIs.YLABEL = 'Frequencies (Hz)';
UDIs.TITLE = {spikeGroups{1};spikeGroups{2}};
UDIs.COLORBARTITLE = 'Coherence';
UDIs.SAVENAME = [saveDirectory 'coherence_matrix_' session '_2groups'];

UDIs = plotAnalyses(dataCell,freqs(freqIdx),UDIs,'matrix');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Histogram of coherence of interacting behaviors by frequency %%%%%%%%%%

%%% 2 Groups %%%
[r1,~] = find((behaviorMatAll(:,:,spikeGroupIdx{1}) == 5) .* (coherenceAll(freqIdx,:,spikeGroupIdx{1}) > 0.4));
[r2,~] = find((behaviorMatAll(:,:,spikeGroupIdx{2}) == 5) .* (coherenceAll(freqIdx,:,spikeGroupIdx{2}) > 0.4));
dataCell{1} = r1; dataCell{2} = r2; 
UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.XAXIS = freqs(freqIdx);
UDIs.TITLE = {'Coherence Distribution'};
UDIs.SAVENAME = [saveDirectory 'coherence_histogram_' session '_2groups'];

UDIs = plotAnalyses(dataCell,freqs(freqIdx),UDIs,'histogramPlot');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Counts of certain behaviors by coherence level %%%%%%%%%%%%%%%%%%%%%%%%

%%% 2 Groups %%%

% NOTE: This plot has a special log colorbar scale!
tk = logspace(0,0.79,25);
cmap = jet(24);
ctable = [tk(1:end-1)' cmap*255 tk(2:end)' cmap*255];
save mycol.cpt ctable -ascii;

figure('Position',[100 100 560*2 500]);
colormap(jet)
handleAX(1) = subplot(121);
[a,aI] = sortrows(sum(histCounts(:,:,spikeGroupIdx{1}),3),1);
imagesc(log(a))
yLabel = {'M';'NM';'S';'SNOT/SNOB';'ITO/IBO';'ATO/ABO';'LTO/LBO'};
set(gca,'XTick',1:9,'XTickLabel',0.1:0.1:0.9,'YTick',1:7,'YTickLabel',yLabel(aI))
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
cptcmap('mycol', 'mapping', 'direct');
hC = cptcbar(gca, 'mycol', 'eastoutside', true);
title(spikeGroups{1})

handleAX(2) = subplot(122);
%[b,bI] = sortrows(sum(histCounts(:,:,spikeGroupIdx{2}),3),1);
b = sum(histCounts(:,:,spikeGroupIdx{2}),3);
h2 = imagesc(log(b(aI,:)));
set(gca,'XTick',1:9,'XTickLabel',0.1:0.1:0.9,'YTick',1:7,'YTickLabel',yLabel(aI))
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
cptcmap('mycol', 'mapping', 'direct');
hC = cptcbar(gca, 'mycol', 'eastoutside', true);
title(spikeGroups{2})

save2pdf([saveDirectory 'coherence_2Dhistogram_' session '_2groups'])

%handleAX(1).Position = [handleAX(1).Position(1) 0.12 (400/(560*2)) 0.8];
%handleAX(2).Position = [handleAX(2).Position(1) 0.12 (400/(560*2)) 0.8];
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scatter with coherence and performance
% [performanceNOR,~]  = readPerformance(codesFile,5,7,6,5:11,analysisAnimals);
% performanceNOR(performanceNOR == 0) = NaN;
% performanceNOR(performanceNOR == 1) = NaN;
% 
% vals = squeeze((behaviorMatAll(6,:,:) == 5 | behaviorMatAll(6,:,:) == 6 | behaviorMatAll(6,:,:) == 7) .* coherenceAll(9,:,:)); % Coherence for interacting behaviors at 7.8 Hz
% groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};
% spikeGroupColor = {[0 0 0];[0.5 0.5 0.5]};
% 
% figure('Position',[100 100 2000 1000])
% subplot(121)
% hold on
% for iGroup = 1:length(groupIdx)
%     for iAnimal = 1:length(groupIdx{iGroup})
%         
%         temp = vals(:,groupIdx{iGroup}(iAnimal));
%         temp = temp((temp ~= 0) & ~isnan(temp));
%         plot(temp,performanceNOR(groupIdx{iGroup}(iAnimal))*ones(1,length(temp)),'.','MarkerSize',20,'Color',groupColor{iGroup})
%         
%     end
% end
% xlabel('Coherence')
% ylabel('NOR Performance')
% subplot(122)
% hold on
% for iGroup = 1:length(spikeGroupIdx)
%     for iAnimal = 1:length(spikeGroupIdx{iGroup})
%         
%         temp = vals(:,spikeGroupIdx{iGroup}(iAnimal));
%         temp = temp((temp ~= 0) & ~isnan(temp));
%         plot(temp,performanceNOR(spikeGroupIdx{iGroup}(iAnimal))*ones(1,length(temp)),'.','MarkerSize',20,'Color',spikeGroupColor{iGroup})
%         
%     end
% end
% xlabel('Coherence')
% ylabel('NOR Performance')

%
% figure('Position',[100 100 500 1000])
% hold on
% temp = behaviorMatAll(:,:,groupIdx{1}) .* (coherenceAll(freqIdx,:,groupIdx{1}) > 0);
% temp(temp == 0) = [];
% histogram(temp(:))
% temp = behaviorMatAll(:,:,groupIdx{2})  .* (coherenceAll(freqIdx,:,groupIdx{2}) > 0);
% temp(temp == 0) = [];
% histogram(temp(:))
% temp = behaviorMatAll(:,:,groupIdx{3}) .* (coherenceAll(freqIdx,:,groupIdx{3}) > 0);
% temp(temp == 0) = [];
% histogram(temp(:))
% temp = behaviorMatAll(:,:,groupIdx{4})  .* (coherenceAll(freqIdx,:,groupIdx{4}) > 0);
% temp(temp == 0) = [];
% histogram(temp(:))

soundAlarm