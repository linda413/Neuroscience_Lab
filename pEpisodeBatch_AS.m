%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Batch P-Episode for Gurkoff Lab 
%
% Amber Schedlbauer (2018)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% User Defined Inputs

% Path for folder where the data is located (i.e. all the files to input
% into pEpisode)
directoryPath = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Converted Files/';

% Name of directory you would specifically like to analyze
dataDirectory = 'A4_novel1_070116';

% Enter a vector of frequencies (log-space) over which you would like to
% calculate pEpisode (entire spectrum)
% Example: frequencies = logspace(log10(3),log10(54),24);
backgroundFrequencies = logspace(log10(2),log10(54),30);

% Enter the frequency band of particular interest to be output to the .csv
% and figure files
frequenciesOfInterest = [6 10];

% Enter the time window in seconds that will be used for the background
% spectrum
% Example: time window 0 min - 5 min -> backgroundWindow = [0 300];
backgroundWindow = [0 300];

% Enter the time window in seconds that you are curious to examine pEpisode
% Example: time window 1 min - 3 min -> windowOfInterest = [60 180]
windowOfInterest = [0 120];

% Path for folder where the results will be saved
saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/' dataDirectory '/'];

%% Pepisode Batch 

% Grabs the files in the data directory
allFiles = rdir([directoryPath dataDirectory '/*.mat']); % FYI rdir is a non-Matlab function

for iFile = 1:length(allFiles)
    
    LFP = load(allFiles(iFile).name);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% INSERT SPIKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insert spike information into LFP structure as bad_intervals
    if isfield(LFP,'spike_intervals') % Check to see if this field exists
        LFP.bad_intervals = [LFP.bad_intervals; LFP.spike_intervals];
    else
        continue % If the field does not exist, then the channel was proabably determined to be artifact.
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract names for CSV-file
    [~,name,~] = fileparts(allFiles(iFile).name);
    name = strrep(name,' ','_'); %Gets rid of any spaces and replaces with underscore
    colNames{iFile} = strrep(name,'-','_'); %Gets rid of any dashes and replaces with underscore
    
    % Calculate Pepisode
    Pepisode(:,iFile) = pEpisode_AS(LFP,backgroundFrequencies,frequenciesOfInterest,backgroundWindow,windowOfInterest,saveDirectory);   
    
    clear LFP removeInterval name
    
end

%% Save outputs to CSV file

% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir') 
    mkdir(saveDirectory)
end

T = array2table([backgroundFrequencies' Pepisode],'VariableNames',[{'Frequencies'} colNames]);
writetable(T,[saveDirectory 'Pepisode_allFiles.csv'])

