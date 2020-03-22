%% Manuscript 2 - Figure 1
%
% EEG traces during NOR 1 in LV
% 1) SHAM - A39
% 2) Low spike - A74 STIM
% 3) High spike - A15 STIM
%
% Need both zoomed in and zoomed out traces

%% Figure 1a

% Simplified timeline of experimental paradigm created in Illustrator

%% Figure 1b

clc; clear; close all

% Load data from 3 animals
sham = load('/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/A39_novel1_071417/A39_novel1_071417_07-Mar-2018_CSC3.mat');
lowSpike = load('/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/A74_novel1_090718/A74_novel1_090718_13-Nov-2018_CSC3.mat');
highSpike = load('/Users/amberschedlbauer/Desktop/Gurkoff Lab/Converted Files/NOR/A15_novel1_091616/A15_novel1_091616_07-Mar-2018_CSC12.mat');

% Raw signal that has been filtered for 60 Hz line noise
d_60 = designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,'DesignMethod','butter','SampleRate',1000);
valuesSham = filtfilt(d_60,sham.values);
valuesLowSpike = filtfilt(d_60,lowSpike.values);
valuesHighSpike = filtfilt(d_60,highSpike.values);

figure('Position',[100 100 1000 800])
subplot(311)
plot(valuesSham(20000:70000),'k')
ylim([min(valuesSham(20000:70000)) max(valuesSham(20000:70000))])
set(gca,'visible','off')
subplot(312)
plot(valuesLowSpike(20000:70000),'k')
ylim([min(valuesLowSpike(20000:70000)) max(valuesLowSpike(20000:70000))])
set(gca,'visible','off')
subplot(313)
plot(valuesHighSpike(20000:70000),'k')
ylim([min(valuesHighSpike(20000:70000)) max(valuesHighSpike(20000:70000))])
set(gca,'visible','off')
save2pdf(['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/EEGtraces' ])

figure('Position',[100 100 1000 800])
subplot(311)
plot(valuesSham(23000:28000),'k')
%ylim([min(valuesSham(23000:28000)) max(valuesSham(23000:28000))])
ylim([min(valuesSham(20000:70000)) max(valuesSham(20000:70000))])
set(gca,'visible','off')
subplot(312)
plot(valuesLowSpike(23000:28000),'k')
%ylim([min(valuesLowSpike(23000:28000)) max(valuesLowSpike(23000:28000))])
ylim([min(valuesLowSpike(20000:70000)) max(valuesLowSpike(20000:70000))])
set(gca,'visible','off')
subplot(313)
plot(valuesHighSpike(23000:28000),'k')
%ylim([min(valuesHighSpike(23000:28000)) max(valuesHighSpike(23000:28000))])
ylim([min(valuesHighSpike(20000:70000)) max(valuesHighSpike(20000:70000))])
set(gca,'visible','off')
save2pdf(['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/EEGtracesZoom' ])

%% Figure 1c

%plotNORbehavior