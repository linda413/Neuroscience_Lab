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
directoryPath = '/shared/Katelynn/20495/Converted Files/Behavior/NOR/';

% Enter either:
%   1) Excel file containing the animal information
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx';
%   2) Mat file containing the "fileInfo" variable (i.e. files2analyze_date.mat)
%      Ex: codesFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/files2analyze_4_2_18.mat';
codesFile = '/shared/Katelynn/20495/20495_Channels.xlsx';

% Task to analyze (enter only 1)
session = 'NOR2';

% Group to analyze
groupInput = {'Sham';'Pilo'};

% Treatment to analyze
treatmentInput = {'None';'Stim'};

% Region to analyze (enter only 1)
electrode = {'dHPC1'};

% Enter a vector of frequencies (log-space) over which you would like to
% calculate power (entire spectrum)
% Example: frequencies = logspace(log10(3),log10(54),24);
frequencies = logspace(log10(2),log10(30),30);

% Enter the frequency band of particular interest
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that you are curious to examine
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [10 25];

% Enter a "1" or "0" if you want to use baseline normalization
% If so, modify the baselineNormalization function to suit your needs.
% ALI: USE 0 for Baseine files and 1 for actual Barnes sessions
useBaselineNormalization = 0;

% Path for folder where the results will be saved
saveDirectory = ['/shared/Katelynn/20495/Results/Power/' savedate '/' session '/' electrode{1} '/'];

% Calculate power for specific behaviors
addSpecificBehaviors = 0;
videoAnalysisFile = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx';
%behavior = {'NM';'S'};
%behavior = {'M'};
%behavior = {'ITO';'IBO'};

%% Identify files to analyze and calculate power

% Grabs the animals of interest
[~,~,ext] = fileparts(codesFile);
if strcmp(ext,'.xlsx')
    %fileInfo = readAnimalCodes(excelFile,animalsInput,groupInput,treatmentInput,electrodeInput,NORinput,saveDir)
    fileInfo = readAnimalCodes(codesFile,[],groupInput,treatmentInput,electrode,0,[]);
elseif strcmp(ext,'.mat')
    load(codesFile)
end
analysisAnimals = fileInfo(:,1);

% Indices of animals belonging to a treatment group
groups = {'SHAM','PILO','STIM'};
groupIdx{1} = find(strcmp(fileInfo(:,2),'SHAM'));
groupIdx{2} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'NO'));
groupIdx{3} = find(strcmp(fileInfo(:,2),'PILO') & strcmp(fileInfo(:,3),'STIM'));
groupLabel = cell(1,length(analysisAnimals));
groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; 

% Read bad channels
channels = readBadChannels(fileInfo,codesFile,2,session);
sheet = 2;

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
    %  d_60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',sRate);
    %  LFP.values = filtfilt(d_60,LFP.values);
        %fvtool(d_60) %Look at filter reponse
        
        %%%%%%%%%%%%%%%%%%%%%% POWER CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
         OUT(iAnimal) = power_AS(LFP,frequencies,windowOfInterest,0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Save calculated information for all animals
        powerAll(:,:,iAnimal) = OUT(iAnimal).power(:,interestWindowIdx);

        % Conduct baseline normalization on raw power values
        if ~useBaselineNormalization
            
            baselinePower = OUT(iAnimal).power;
            save([directoryPath analysisAnimals{iAnimal} '_baselinePower'],'baselinePower')
            
        elseif useBaselineNormalization
            
            load([directoryPath analysisAnimals{iAnimal} '_baselinePower.mat'])
            
            % Block time window (e.g. baselineWindow = [0 300];)
            % timeStamps = 0:(1/sRate):(size(powerAll(:,:,iAnimal),2)/sRate);
            baselineWindow = [1 290];
            timeVec = logical(timeStamps >= baselineWindow(1) & timeStamps < baselineWindow(2)); % Baseline index
            
            if mean(isnan(baselinePower(1,timeVec))) > 0.5
                warning('Baseline window is more than 50% artifact!')
            end
            
            baselineMean = nanmean(baselinePower(:,timeVec),2);
            
            % Decibel conversion
            powerAllConvert(:,:,iAnimal) = 10*log10(bsxfun(@rdivide,powerAll(:,:,iAnimal),baselineMean));
            
        end
        
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
PSD1 = squeeze(nanmean(powerAll,2));
PSD2 = squeeze(nanmean(powerAllConvert,2));
%powerAll31 = powerAll;

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

%% PLOT AND/OR SAVE

soundAlarm
