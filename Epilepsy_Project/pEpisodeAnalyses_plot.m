%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P-Episode analyses for Gurkoff Lab

% Amber Schedlbauer - 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located (i.e. all the files to input into pEpisode)
%directoryPath = '/Volumes/Gurkoff Lab/Amber/Converted Files/Epileptogenesis/';
directoryPath = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/';
%directoryPath = '/Volumes/SP PHD U3/AmberS/Converted Files/NOR/';

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
% calculate pEpisode (entire spectrum)
% Example: frequencies = logspace(log10(3),log10(54),24);
backgroundFrequencies = logspace(log10(2),log10(30),30);

% Enter the frequency band of particular interest to be output to the .csv and figure files
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that will be used for the background spectrum
% Example: time window 0 min - 5 min -> backgroundWindow = [0 300];
backgroundWindow = [0 300];

% Enter the time window in seconds that you are curious to examine pEpisode
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [0 120];

% Path for folder where the results will be saved
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/' session '/' electrode{1} '/'];

% Calculate Pepidose for specific behaviors
addSpecificBehaviors = 1;
videoAnalysisFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx';
%behavior = {'NM';'S'};
%behavior = {'M'};
behavior = {'ITO';'IBO'};

%% Identify files to analyze and calculate Pepisode

% Grabs the animals of interest
[~,~,ext] = fileparts(codesFile);
if strcmp(ext,'.xlsx')
    %fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
    fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrode,1,[]);
elseif strcmp(ext,'.mat')
    load(codesFile)
end
analysisAnimals = fileInfo(:,1);

% Indices of animals belonging to a group
groups = {'SHAM','PILO','STIM','BURST'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupIdx{4} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'BURST'));
groupLabel = cell(1,length(analysisAnimals));
groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; groupLabel(groupIdx{4}) = {'BURST'};

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
[~,idxF1] = min(abs(backgroundFrequencies - frequenciesOfInterest(1)));
[~,idxF2] = min(abs(backgroundFrequencies - frequenciesOfInterest(2)));

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

% Allocations
PepisodeAll = nan(size(fileInfo,1),length(backgroundFrequencies));
freqTimeWinAll = nan(length(backgroundFrequencies),length(windowOfInterest(1):windowOfInterest(2)-1)*1000,size(fileInfo,1));
powerAll = nan(length(backgroundFrequencies),length(backgroundWindow(1):backgroundWindow(2)-1)*1000,size(fileInfo,1));
behaviorPepisode = nan(size(fileInfo,1),length(backgroundFrequencies));
behaviorPepisodeAlign = nan(length(backgroundFrequencies),3000,size(fileInfo,1));
[numWindows,percentArtifactBackground,percentArtifactInterestWindow] = deal(nan(size(fileInfo,1),1));

%% Loop through animals and removes bad channels and computes p-episode
for iAnimal = 1:size(analysisAnimals)
    
    if ~isempty(channels{iAnimal})
        
        % Load file
        fileName = rdir([directoryPath analysisAnimals{iAnimal} '_' session '_*/' analysisAnimals{iAnimal} '_*' channels{iAnimal} '.mat']);
        if exist(fileName.name,'file') == 2
            LFP = load(fileName.name); disp(fileName.name);
        else
            continue
        end
        
        % Get time/freq info
        sRate = LFP.sFreq;
        timeStamps = round(LFP.timestamps - LFP.timestamps(1),3); % Grabs total signal time with timepoint(1) = 0
        interestWindowIdx = find(timeStamps >= windowOfInterest(1) & timeStamps < windowOfInterest(2));
        
        % Insert spike information into LFP structure as bad_intervals
        if isfield(LFP,'spike_intervals') % Check to see if this field exists
            LFP.bad_intervals = [LFP.bad_intervals; LFP.spike_intervals];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% PEPISODE CALCULATIONS %%%%%%%%%%%%%%%%%%%%%
        pEpisodeDir = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Pepisode/';
        pEpisodeFile = [pEpisodeDir 'Pepisode_' analysisAnimals{iAnimal} '_' session '_' channels{iAnimal} '.mat'];
        if exist(pEpisodeFile,'file') == 2
            load(pEpisodeFile)
        else
            % Filter raw signal for 60 Hz line noise
            d_60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',sRate);
            LFP.values = filtfilt(d_60,LFP.values);
            %fvtool(d_60) %Look at filter reponse
            
            [Pepisode,freqTimeWin,power] = pEpisode_AS(LFP,backgroundFrequencies,frequenciesOfInterest,backgroundWindow,windowOfInterest,[]);
            
            if ~exist(pEpisodeDir,'dir')
                mkdir(pEpisodeDir)
            end
            save(pEpisodeFile(1:end-4),'backgroundFrequencies','frequenciesOfInterest','backgroundWindow','windowOfInterest','Pepisode','freqTimeWin','power')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Finds portion of signal that is artifact/spike
        artifactBackground = isnan(power(1,:));
        percentArtifactBackground(iAnimal) = mean(artifactBackground);
        percentArtifactInterestWindow(iAnimal) = mean(artifactBackground(interestWindowIdx));
        
        % Save calculated information for all animals
        PepisodeAll(iAnimal,:) = Pepisode;
        freqTimeWinAll(:,:,iAnimal) = freqTimeWin;
        powerAll(:,:,iAnimal) = power;
        
        %         % Sanity Checks: Plot individual raw data, power, & Pepisode
        %         figure('Position',[100 100 1500 1000])
        %         ax(1) = subplot(311);
        %         backWindowIdx = timeStamps >= backgroundWindow(1) & timeStamps < backgroundWindow(2);
        %         rawData = LFP.values(backWindowIdx); rawData(artifactBackground) = NaN;
        %         plot(timeStamps(interestWindowIdx),rawData(interestWindowIdx))
        %         title(saveName)
        %         ax(2) = subplot(312);
        %         imagesc(timeStamps(interestWindowIdx),1:length(backgroundFrequencies),squeeze(power(:,interestWindowIdx)))
        %         set(gca,'YTick',1:5:length(backgroundFrequencies),'YTickLabel',backgroundFrequencies(1:5:end))
        %         title('Power')
        %         ax(3) = subplot(313);
        %         imagesc(timeStamps(interestWindowIdx),1:length(backgroundFrequencies),squeeze(freqTimeWin))
        %         set(gca,'YTick',1:5:length(backgroundFrequencies),'YTickLabel',backgroundFrequencies(1:5:end))
        %         title('Pepisode')
        %         linkaxes(ax,'x')
        %         %saveName = [analysisAnimals{iAnimal} ' ' analysisChan ' Raw Power Pepisode'];
        %         %saveas(gcf,[saveDirectory saveName])
        
    else
        continue
    end
    
    % Calculates Pepisode for the time windows that correspond to the
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
                freqTimeWinAlign(:,:,iTime) = freqTimeWin(:,startTime:endTime);
                
            end
            
            % Averages across all time points (multiple interacting windows)
            behaviorPepisode(iAnimal,:) = nanmean(freqTimeWin(:,timeVec),2)*100;
            
            % Adds across all interacting windows
            behaviorPepisodeAlign(:,:,iAnimal) = nansum(freqTimeWinAlign,3);
            
        end
        
    end
    
    clear filename LFP pEpisodeFile fileNames timeWindow freqTimeWin freqTimeWinAlign
    
end

%% Plot

% Colors for treatment groups
groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};
spikeGroupColor = {[0 0 0];[0.5 0.5 0.5]};

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

%%% All Behaviors Pepisode "PSD" Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 4 Groups %%%
dataCell{1} = nanmean(PepisodeAll(groupIdx{1},:)); dataCell{2} = nanmean(PepisodeAll(groupIdx{2},:));
dataCell{3} = nanmean(PepisodeAll(groupIdx{3},:)); dataCell{4} = nanmean(PepisodeAll(groupIdx{4},:));

UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.XAXIS = backgroundFrequencies;
UDIs.XLABEL = {'Frequencies (Hz)'};
UDIs.YLABEL = '% Time';
UDIs.SAVENAME = [saveDirectory 'PSDmean_allBehaviors_4groups'];

UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'PSD');
clear data*
%%%

%%% 2 Groups %%%
dataCell{1} = PepisodeAll(spikeGroupIdx{1},:); dataCell{2} = PepisodeAll(spikeGroupIdx{2},:);

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.XAXIS = backgroundFrequencies;
UDIs.XLABEL = {'Frequencies (Hz)'};
UDIs.YLABEL = '% Time';
UDIs.SAVENAME = [saveDirectory 'PSDmean_allBehaviors_2groups'];

UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'PSDsig');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

%%% All Behaviors Average Pepisode Boxplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avePepisode = nanmean(PepisodeAll(:,idxF1:idxF2),2);
avePepisode(percentArtifactInterestWindow > 0.5) = NaN;

% Gene's Revised Pepisode
%avePepisode = squeeze(sum(freqTimeWinAll(idxF1:idxF2,:,:),1));
%nanIdx = isnan(sum(avePepisode));
%avePepisode = squeeze(nanmean(avePepisode > 0))';
%avePepisode(nanIdx) = NaN;

%%% 4 Groups %%%
dataCell = {avePepisode(groupIdx{1})',avePepisode(groupIdx{2})',avePepisode(groupIdx{3})',avePepisode(groupIdx{4})'};

UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.YLABEL = ['Average % Time (' num2str(frequenciesOfInterest(1)) '-' num2str(frequenciesOfInterest(2)) 'Hz)'];

[pval,~,~] = anova1(avePepisode,groupLabel,'off'); % Unbalanced 1-way ANOVA for significance
UDIs.TITLE = {['1-Way ANOVA: p = ' num2str(pval)]};

UDIs.SAVENAME = [saveDirectory 'pepisode_box_allBehaviors_4groups'];

UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'boxPlotScatter');
clear data*
%%%

%%% 2 Groups %%%
dataCell = {avePepisode(spikeGroupIdx{1})',avePepisode(spikeGroupIdx{2})'};

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.YLABEL = ['Average % Time (' num2str(frequenciesOfInterest(1)) '-' num2str(frequenciesOfInterest(2)) 'Hz)'];

[~,pval] = ttest2(dataCell{1},dataCell{2}); % Non-paired ttest for significance
UDIs.TITLE = {['2-Sample Ttest: p = ' num2str(pval)]};

UDIs.SAVENAME = [saveDirectory 'pepisode_box_allBehaviors_2groups'];

UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'boxPlotScatter');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% All Behaviors Pepisode Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 4 Groups %%%
dataCell{1} = nanmean(freqTimeWinAll(:,:,groupIdx{1}),3); dataCell{2} = nanmean(freqTimeWinAll(:,:,groupIdx{2}),3);
dataCell{3} = nanmean(freqTimeWinAll(:,:,groupIdx{3}),3); dataCell{4} = nanmean(freqTimeWinAll(:,:,groupIdx{4}),3);

UDIs.groups = groups;
UDIs.groupColor = groupColor;
UDIs.XAXIS = timeStamps(interestWindowIdx);
UDIs.YAXIS = backgroundFrequencies;
UDIs.COLORMAP = flipud(gray);
UDIs.XLABEL = {'','','Time (s)',''};
UDIs.YLABEL = 'Frequencies (Hz)';
UDIs.TITLE = groups;
UDIs.COLORBARTITLE = 'Percent';
UDIs.SAVENAME = [saveDirectory 'pepisode_matrix_allBehaviors_4groups'];

UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'matrix');
clear data*
%%%

%%% 2 Groups %%%%
dataCell{1} = nanmean(freqTimeWinAll(:,:,spikeGroupIdx{1}),3); dataCell{2} = nanmean(freqTimeWinAll(:,:,spikeGroupIdx{2}),3);

UDIs.groups = spikeGroups;
UDIs.groupColor = spikeGroupColor;
UDIs.XAXIS = timeStamps(interestWindowIdx);
UDIs.YAXIS = backgroundFrequencies;
UDIs.COLORMAP = flipud(gray);
UDIs.XLABEL = {'Time (s)','Time (s)'};
UDIs.YLABEL = 'Frequencies (Hz)';
UDIs.TITLE = spikeGroups;
UDIs.COLORBARTITLE = 'Percent';
UDIs.SAVENAME = [saveDirectory 'pepisode_matrix_allBehaviors_2groups'];

UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'matrix');
clear data*
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Specific Behaviors
if addSpecificBehaviors
    
    %%% All Behaviors Pepisode "PSD" Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 4 Groups %%%
    dataCell{1} = nanmean(behaviorPepisode(groupIdx{1},:)); dataCell{2} = nanmean(behaviorPepisode(groupIdx{2},:));
    dataCell{3} = nanmean(behaviorPepisode(groupIdx{3},:)); dataCell{4} = nanmean(behaviorPepisode(groupIdx{4},:));
    
    UDIs.groups = groups;
    UDIs.groupColor = groupColor;
    UDIs.XAXIS = backgroundFrequencies;
    UDIs.XLABEL = {'Frequencies (Hz)'};
    UDIs.YLABEL = '% Time';
    UDIs.SAVENAME = [saveDirectory 'PSDmean_4groups'];
    
    UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'PSD');
    clear data*
    %%%
    
    %%% 2 Groups %%%
    dataCell{1} = behaviorPepisode(spikeGroupIdx{1},:); dataCell{2} = behaviorPepisode(spikeGroupIdx{2},:);
    
    UDIs.groups = spikeGroups;
    UDIs.groupColor = spikeGroupColor;
    UDIs.XAXIS = backgroundFrequencies;
    UDIs.XLABEL = {'Frequencies (Hz)'};
    UDIs.YLABEL = '% Time';
    UDIs.SAVENAME = [saveDirectory 'PSDmean_2groups'];
    
    UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'PSDsig');
    clear data*
    %%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Specific Behaviors Average Pepisode Boxplot %%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 4 Groups %%%
    dataCell = {behaviorPepisode(groupIdx{1})',behaviorPepisode(groupIdx{2})',behaviorPepisode(groupIdx{3})',behaviorPepisode(groupIdx{4})'};
    
    UDIs.groups = groups;
    UDIs.groupColor = groupColor;
    UDIs.YLABEL = ['Average % Time (' num2str(frequenciesOfInterest(1)) '-' num2str(frequenciesOfInterest(2)) 'Hz)'];
    
    [pval,~,~] = anova1(behaviorPepisode,groupLabel,'off'); % Unbalanced 1-way ANOVA for significance
    UDIs.TITLE = {['1-Way ANOVA: p = ' num2str(pval)]};
    
    UDIs.SAVENAME = [saveDirectory 'pepisode_box_4groups'];
    
    UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'boxPlotScatter');
    clear data*
    %%%
    
    %%% 2 Groups %%%
    dataCell = {behaviorPepisode(spikeGroupIdx{1})',behaviorPepisode(spikeGroupIdx{2})'};
    
    UDIs.groups = spikeGroups;
    UDIs.groupColor = spikeGroupColor;
    UDIs.YLABEL = ['Average % Time (' num2str(frequenciesOfInterest(1)) '-' num2str(frequenciesOfInterest(2)) 'Hz)'];
    
    [~,pval] = ttest2(dataCell{1},dataCell{2}); % Non-paired ttest for significance
    UDIs.TITLE = {['2-Sample Ttest: p = ' num2str(pval)]};
    
    UDIs.SAVENAME = [saveDirectory 'pepisode_box_2groups'];
    
    UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'boxPlotScatter');
    clear data*
    %%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Specific Behaviors Pepisode Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 4 Groups %%%
    dataCell{1} = squeeze(nansum(behaviorPepisodeAlign(:,:,groupIdx{1}),3))/nansum(numWindows(groupIdx{1}));
    dataCell{2} = squeeze(nansum(behaviorPepisodeAlign(:,:,groupIdx{2}),3))/nansum(numWindows(groupIdx{2}));
    dataCell{3} = squeeze(nansum(behaviorPepisodeAlign(:,:,groupIdx{3}),3))/nansum(numWindows(groupIdx{3}));
    dataCell{4} = squeeze(nansum(behaviorPepisodeAlign(:,:,groupIdx{4}),3))/nansum(numWindows(groupIdx{4}));
    
    UDIs.groups = groups;
    UDIs.groupColor = groupColor;
    UDIs.XAXIS = -0.999:0.001:2;
    UDIs.YAXIS = backgroundFrequencies;
    UDIs.COLORMAP = flipud(gray);
    UDIs.XLABEL = {'','','Time (s)',''};
    UDIs.YLABEL = 'Frequencies (Hz)';
    UDIs.TITLE = groups;
    UDIs.COLORBARTITLE = 'Percent';
    UDIs.SAVENAME = [saveDirectory 'pepisode_matrix_4groups'];
    
    UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'matrix');
    clear data*
    %%%
    
    %%% 2 Groups %%%%
    dataCell{1} = squeeze(nansum(behaviorPepisodeAlign(:,:,spikeGroupIdx{1}),3))/nansum(numWindows(spikeGroupIdx{1}));
    dataCell{2} = squeeze(nansum(behaviorPepisodeAlign(:,:,spikeGroupIdx{2}),3))/nansum(numWindows(spikeGroupIdx{2}));
    
    UDIs.groups = spikeGroups;
    UDIs.groupColor = spikeGroupColor;
    UDIs.XAXIS = -0.999:0.001:1.999;
    UDIs.YAXIS = backgroundFrequencies;
    UDIs.COLORMAP = flipud(gray);
    UDIs.XLABEL = {'Time (s)','Time (s)'};
    UDIs.YLABEL = 'Frequencies (Hz)';
    UDIs.TITLE = spikeGroups;
    UDIs.COLORBARTITLE = 'Percent';
    UDIs.SAVENAME = [saveDirectory 'pepisode_matrix_2groups'];
    
    UDIs = plotAnalyses(dataCell,backgroundFrequencies,UDIs,'matrix');
    clear data*
    %%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

soundAlarm