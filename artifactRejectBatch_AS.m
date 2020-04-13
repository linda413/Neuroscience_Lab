%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to call DAVIS_MultiTrace_Manual_Artifact_Reject and input files
% rather than use the GUI -> makes it a little faster
%
% Amber Schedlbauer
% 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% UDIs

% Enter the directory where the data are located
%mainDir = '/Volumes/SP PHD U3/AmberS/Converted Files/NOR/';
%mainDir = '/home/amschedl/Desktop/Amber/Converted Files/Epileptogenesis/';
mainDir = '/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/';
mainDir = '/Users/caiqi/Desktop/Neuroscience_Lab/'

% Enter the recording session you would like to artifact reject
session = 'BMD';

% Enter the prefix used for animal identification
animalPrefix = 'HM'; % Change the "A" if your animals are named differently!
channels = [1,4,5,8,12,13,16] ;
% Enter the number of seconds you would like to view per window
numSec = 1450;

%% Loop through animals

directories = dir([mainDir animalPrefix '*' session '*']);

for i = 1:length(directories)
    
    
    % Identifies which channels will be viusalized and artifact rejected
    newDir = directories(i).name;
    

    allFiles = {};
    filesPath = {};
    for j = 1: length(channels)
        channel = channels(j)
        getChannel = append('*CSC',string(channel),'.mat')
        file = dir([mainDir newDir '/' convertStringsToChars(getChannel)])
        allFiles = [allFiles file.name];
        
 

    % Call the function to simultaneously view all channels for artifact rejection
    end
    filePath = [mainDir newDir '/' ]
    filesPath = filePath

    disp(allFiles)
    DAVIS_MultiTrace_Manual_Artifact_Reject(allFiles,filesPath,numSec,0)

    % Pauses the loop until the figure window is closed
    waitfor(gcf)   
    
    
    
end