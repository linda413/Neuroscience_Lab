%% seeAllChannels

clc; clear; close all;

%% UDIs

% Enter the path of where all the data is located
%mainDir = '/home/amschedl/Desktop/NOR/ConvertedFiles/';
%mainDir = '/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Converted Files/';
mainDir = '/Volumes/SP PHD U3/AmberS/Converted Files/';

% Enter the animal you want to analyze
animalNum = 13;

% Enter task you want to analyze
task = 'novel1';

%% Grab filenames

% Load output from readNORcodes.m
%load('/home/amschedl/Desktop/NOR/Code/Results/files2analyze_3_1_18.mat')
load('/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/files2analyze_4_2_18.mat')

files = rdir([mainDir fileInfo{animalNum,1} '_' task '_*/*.mat']);
allFiles = {files.name}';

figure('Position',[100 100 1500 750])

for iPlot = 1:length(allFiles)
    
    LFP = load(allFiles{iPlot});
    
    ts = LFP.timestamps - LFP.timestamps(1);
    %ax(iPlot) = subplot(4,1,iPlot)
    plot(ts,LFP.values,'k')
    title(allFiles{iPlot})
    hold on
    pause
    
end

%close(gcf)
%linkaxes([ax(1) ax(2) ax(3) ax(4)],'xy')
