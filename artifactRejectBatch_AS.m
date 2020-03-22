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

% Enter the recording session you would like to artifact reject
session = 'novel2';

% Enter the prefix used for animal identification
animalPrefix = 'A'; % Change the "A" if your animals are named differently!

% Enter the number of seconds you would like to view per window
numSec = 450;

%% Loop through animals

directories = dir([mainDir animalPrefix '*' session '*']);

for i = 1:length(directories)
    
    % Identifies which channels will be viusalized and artifact rejected
    newDir = directories(i).name;
    allFiles = rdir([mainDir newDir '/']);
    allLFPs = {allFiles(:).name}';
    disp(allLFPs) 

    % Call the function to simultaneously view all channels for artifact rejection
    DAVIS_MultiTrace_Manual_Artifact_Reject(allLFPs,numSec,1)
    
    % Pauses the loop until the figure window is closed
    waitfor(gcf)
    
end
 