%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Batch Conversion of NSC to MAT files for Gurkoff Lab
%
% Amber Schedlbauer
% 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located (i.e. all the files to input)
% If you would like to analyze multiple folders at a time, please enter the
% directory where those folders are located.  
dataDirectory = '/home/amschedl/Desktop/NOR/EEG Data/Epileptogenesis/';

% Path for folder where the results will be saved
%mainSaveDir = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Code/Results/';
saveDirectory = '/Volumes/SP PHD U3/AmberS/Converted Files/';

% Enter either:
%   1) channels you would like to analyze (Ex: 3 for CSC3.ncs)
%      Ex: channelsToAnalyze = 3;
%   2) the path and name of the MAT file that contains the varibale
%      "fileInfo" (i.e. files2analyze_date.mat)
%      Ex: channelsToAnalyze = '/Users/amberschedlbauer/files2analyze_4_2_18.mat';
channelsToAnalyze = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/files2analyze_4_2_18.mat';

% Prefix for the animal names (animals in this lab are usually labeled prefix + number, e.g. A54)
animalPrefix = 'A'; 

% The data are collected at 3200 Hz; this is plenty for theta analyses, so
% we downsample to 1000 Hz to make analysis faster.
newSamplingFrequency = 1000;

%% Batch Conversion

% Find the names of the folders in the main data directory
dirName = dir([dataDirectory animalPrefix '*']);
dirName = {dirName.name};

% Loop through the folders
for iDir = 1:length(dirName)
   
    % Skips this folder if it already exists (saves on computation time)
    tempSaveDirectory = [saveDirectory dirName{iDir} '/'];
    if exist(tempSaveDirectory,'dir') == 7
        continue
    end
    
    % Coverts specified channels above to file names OR
    % Identifies only those files that should have signal
    if ischar(channelsToAnalyze)
        load(channelsToAnalyze) % Loads the fileInfo file
        animalName = strtok(dirName{iDir},'_');
        animalIdx = strcmp(fileInfo(:,1),animalName);
        filesToAnalyze = fileInfo(animalIdx,4:end);
        filesToAnalyze(cellfun(@isempty,filesToAnalyze)) = []; % remove empty channels
        filesToAnalyze = cellfun(@(x) [x '.ncs'],filesToAnalyze,'UniformOutput',false); %add on file extension
    else
        for iChan = 1:length(channelsToAnalyze)
            % If you want to grab a different CSC file, you could add in
            % extra code below. Ex: filesToAnalyze{iChan} = ['CSC' num2str(channelsToAnalyze(iChan)) '_001.ncs']
            filesToAnalyze{iChan} = ['CSC' num2str(channelsToAnalyze(iChan)) '.ncs'];
        end
    end
    
    %%% Checks file size - Removes if too small %%%
    % Comment this out if it doesn't matter to you
    removeIdx = [];
    for iChan = 1:length(filesToAnalyze)
        tempFile = dir([dataDirectory dirName{iDir} '/' filesToAnalyze{iChan}]);
        if tempFile.bytes < 200000
            removeIdx = [removeIdx iChan];
        end
    end
    filesToAnalyze(removeIdx) = [];
    filesToAnalyze = {filesToAnalyze.name};
    
    %%%%%%%%%%%%%% Main function that converts the NCS files %%%%%%%%%%%%%%
    convertNCS2MAT_AS([dataDirectory dirName{iDir} '/'],saveDirectory,filesToAnalyze,newSamplingFrequency,[]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end