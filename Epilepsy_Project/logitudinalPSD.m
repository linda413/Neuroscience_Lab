clc; clear; close all;

groupColor = {[0 0 1];[1 0 0];[0 1 0];[1 0 1]};

saveDirectory = ['/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Results/' savedate '/'];
% Makes saveDirectory if it doesn't exist
if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory)
end

%% Baseline normalization to power of prepilostim 01 02 03 %%%%%%%%%%%%%%%

% load('/Users/amberschedlbauer/Desktop/Gurkoff Lab/PowerValuespre_raw.mat')
% load('/Users/amberschedlbauer/Desktop/Gurkoff Lab/PowerValues4_raw.mat')
% load('/Users/amberschedlbauer/Desktop/Gurkoff Lab/PowerValues16_raw.mat')
% load('/Users/amberschedlbauer/Desktop/Gurkoff Lab/PowerValues31_raw.mat')
%
% powerPre = cat(2,powerAll1,powerAll2,powerAll3);
%
% baselineMean = squeeze(nanmean(powerPre,2));
% baselineSTD = squeeze(nanstd(powerPre,[],2));
%
% for iAnimal = 1:50
%
%     temp = powerAll4(:,:,iAnimal);
%     powerZ = bsxfun(@minus,temp,baselineMean(:,iAnimal));
%     powerZ = bsxfun(@rdivide,powerZ,baselineSTD(:,iAnimal));
%     psdZ4(:,iAnimal) = nanmean(powerZ,2);
%     clear temp powerZ
%
%     temp = powerAll16(:,:,iAnimal);
%     powerZ = bsxfun(@minus,temp,baselineMean(:,iAnimal));
%     powerZ = bsxfun(@rdivide,powerZ,baselineSTD(:,iAnimal));
%     psdZ16(:,iAnimal) = nanmean(powerZ,2);
%     clear temp powerZ
%
%     temp = powerAll31(:,:,iAnimal);
%     powerZ = bsxfun(@minus,temp,baselineMean(:,iAnimal));
%     powerZ = bsxfun(@rdivide,powerZ,baselineSTD(:,iAnimal));
%     psdZ31(:,iAnimal) = nanmean(powerZ,2);
%     clear temp powerZ
%
% end
%
% PSD = cat(3,psdZ4,psdZ16,psdZ31);
%
% PSDall{1} = squeeze(nanmean(PSD(:,groupIdx{1},:),2));
% PSDall{2}  = squeeze(nanmean(PSD(:,groupIdx{2},:),2));
% PSDall{3}  = squeeze(nanmean(PSD(:,groupIdx{3},:),2));
% PSDall{4}  = squeeze(nanmean(PSD(:,groupIdx{4},:),2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Other

load('/Users/amberschedlbauer/Desktop/Gurkoff Lab/longitudinalPSD/PSDvalues_dB.mat')

PSDpre = nanmean(cat(3,PSD1,PSD2,PSD3),3);

PSD = cat(3,PSDpre,PSD4,PSD16,PSD31);

PSDall{1} = squeeze(nanmean(PSD(:,groupIdx{1},:),2));
PSDall{2}  = squeeze(nanmean(PSD(:,groupIdx{2},:),2));
PSDall{3}  = squeeze(nanmean(PSD(:,groupIdx{3},:),2));
PSDall{4}  = squeeze(nanmean(PSD(:,groupIdx{4},:),2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot

% Index of range of frequencies of interest
[~,idxF1] = min(abs(frequencies - 6));
[~,idxF2] = min(abs(frequencies - 10));

%%% All mean lines on one graph
figure('Position',[100 100 2000 1000])
ylabel('Raw Power')
xlabel('Frequencies (Hz)')

subplot(2,4,1)
hold on
plot(frequencies,PSDall{1}(:,1),'-','Color',groupColor{1},'LineWidth',2);
plot(frequencies,PSDall{1}(:,2),'--','Color',groupColor{1},'LineWidth',2)
plot(frequencies,PSDall{1}(:,3),':','Color',groupColor{1},'LineWidth',2)
plot(frequencies,PSDall{1}(:,4),'-.','Color',groupColor{1},'LineWidth',2)
legend({'PrePilo','PostPilo4','PostPilo16','PostPilo31'})

subplot(2,4,2)
hold on
plot(frequencies,PSDall{2}(:,1),'-','Color',groupColor{2},'LineWidth',2);
plot(frequencies,PSDall{2}(:,2),'--','Color',groupColor{2},'LineWidth',2)
plot(frequencies,PSDall{2}(:,3),':','Color',groupColor{2},'LineWidth',2)
plot(frequencies,PSDall{2}(:,4),'-.','Color',groupColor{2},'LineWidth',2)
legend({'PrePilo','PostPilo4','PostPilo16','PostPilo31'})

subplot(2,4,5)
hold on
plot(frequencies,PSDall{3}(:,1),'-','Color',groupColor{3},'LineWidth',2);
plot(frequencies,PSDall{3}(:,2),'--','Color',groupColor{3},'LineWidth',2)
plot(frequencies,PSDall{3}(:,3),':','Color',groupColor{3},'LineWidth',2)
plot(frequencies,PSDall{3}(:,4),'-.','Color',groupColor{3},'LineWidth',2)
legend({'PrePilo','PostPilo4','PostPilo16','PostPilo31'})

subplot(2,4,6)
hold on
plot(frequencies,PSDall{4}(:,1),'-','Color',groupColor{4},'LineWidth',2);
plot(frequencies,PSDall{4}(:,2),'--','Color',groupColor{4},'LineWidth',2)
plot(frequencies,PSDall{4}(:,3),':','Color',groupColor{4},'LineWidth',2)
plot(frequencies,PSDall{4}(:,4),'-.','Color',groupColor{4},'LineWidth',2)
legend({'PrePilo','PostPilo4','PostPilo16','PostPilo31'})

subplot(2,4,[3 4 7 8])
hold on
for iGroup = 1:length(PSDall)
    
    handleLine(iGroup) = plot(frequencies,PSDall{iGroup}(:,1),'-','Color',groupColor{iGroup},'LineWidth',2);
    plot(frequencies,PSDall{iGroup}(:,2),'--','Color',groupColor{iGroup},'LineWidth',2)
    plot(frequencies,PSDall{iGroup}(:,3),':','Color',groupColor{iGroup},'LineWidth',2)
    plot(frequencies,PSDall{iGroup}(:,4),'-.','Color',groupColor{iGroup},'LineWidth',2)
    
end

legend(handleLine,groups)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Divided by group
figure('Position',[100 100 2000 500])
titles = groups;

for iPlot = 1:4
    
    subplot(1,4,iPlot)
    hold on
    plot(frequencies,PSDall{iPlot}(:,1),'-','Color',groupColor{iPlot},'LineWidth',2);
    plot(frequencies,PSDall{iPlot}(:,2),'--','Color',groupColor{iPlot},'LineWidth',2)
    plot(frequencies,PSDall{iPlot}(:,3),':','Color',groupColor{iPlot},'LineWidth',2)
    plot(frequencies,PSDall{iPlot}(:,4),'-.','Color',groupColor{iPlot},'LineWidth',2)
    
    if iPlot == 1
        ylabel('Power (dB)')
    end
    
    title(titles{iPlot})
    legend({'PrePilo','PostPilo4','PostPilo16','PostPilo31'})
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    
end

save2pdf([saveDirectory 'longitudinalPSD_group'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Divided by session
figure('Position',[100 100 2000 500])
titles = {'PrePilo','PostPilo4','PostPilo16','PostPilo31'};
lineStyle = {'-','--',':','-.'};

for iPlot = 1:4
    
    subplot(1,4,iPlot)
    hold on
    plot(frequencies,PSDall{1}(:,iPlot),'Color',groupColor{1},'LineWidth',2,'LineStyle',lineStyle{iPlot})
    plot(frequencies,PSDall{2}(:,iPlot),'Color',groupColor{2},'LineWidth',2,'LineStyle',lineStyle{iPlot})
    plot(frequencies,PSDall{3}(:,iPlot),'Color',groupColor{3},'LineWidth',2,'LineStyle',lineStyle{iPlot})
    plot(frequencies,PSDall{4}(:,iPlot),'Color',groupColor{4},'LineWidth',2,'LineStyle',lineStyle{iPlot})
    %ylim([0 0.035])
    if iPlot == 1
        ylabel('Power (dB)')
    end
    title(titles{iPlot})
    legend(groups)
    set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
    
end

save2pdf([saveDirectory 'longitudinalPSD_session'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% All lines over 16 plots
ctr = 1;
figure('Position',[100 100 2000 1000])
for iPlot = 1:size(PSD,3)
    
    for iGroup = 1:4
        
        subplot(size(PSD,3),4,ctr)
        plot(frequencies,squeeze(PSD(:,groupIdx{iGroup},iPlot)),'LineWidth',0.5)
        hold on
        plot(frequencies,PSDall{iGroup}(:,iPlot),'Color',groupColor{iGroup},'LineWidth',4)
        %ylim([-10 0])
        
        ctr = ctr + 1;
        
    end
end

%%% Single value for theta band power

theta = squeeze(nanmean(PSD(idxF1:idxF2,:,:),1));
%theta = squeeze(max(PSD(idxF1:idxF2,:,:)));

figure('Position',[100 100 1000 1000])
hold on

for iGroup = 1:length(groups)
    
    plot(1:size(PSD,3),theta(groupIdx{iGroup},1:size(PSD,3)),'Color',[groupColor{iGroup} 0.25],'LineWidth',0.5)
    handleLine(iGroup) = plot(1:size(PSD,3),nanmean(theta(groupIdx{iGroup},1:size(PSD,3)),1),'Color',groupColor{iGroup},'LineWidth',4);
    
end

set(gca,'XTick',1:size(PSD,3),'XTickLabel',{'Pre','Post4','Post16','Post31'})
ylabel('Single Theta Value')
set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
legend(handleLine,groups)
title('rmANOVA: not significant')

% Repeated Measures ANOVA
fileInfo = readAnimalCodes('/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Animal Codes.xlsx',[],{'SHAM';'PILO'},{'NO'; 'STIM';'BURST'},{'LV'},1,[]);
analysisAnimals = fileInfo(:,1);
groupLabel = cell(length(analysisAnimals),1);
groupLabel(groupIdx{1}) = {'SHAM'}; groupLabel(groupIdx{2}) = {'PILO'}; groupLabel(groupIdx{3}) = {'STIM'}; groupLabel(groupIdx{4}) = {'BURST'};
t = table(groupLabel,theta(:,1),theta(:,2),theta(:,3),theta(:,4),'VariableNames',{'group','meas1','meas2','meas3','meas4'});
Meas = table([1 2 3 4]','VariableNames',{'Session'});
rm = fitrm(t,'meas1-meas4~group','WithinDesign',Meas);
ranovatbl = ranova(rm)
%%% There is no group by session interaction.


save2pdf([saveDirectory 'longitudinalTheta'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Percent change
%
% theta = squeeze(nanmean(PSD(idxF1:idxF2,:,:),1));
% percentChange = bsxfun(@minus,theta,theta(:,1));
% percentChange = bsxfun(@rdivide,percentChange,theta(:,1));
%
% figure('Position',[100 100 1000 1000])
% hold on
%
% for iGroup = 1:length(groups)
%
%     plot(1:4,percentChange(groupIdx{iGroup},1:4),'Color',groupColor{iGroup},'LineWidth',1)
%     %plot(1:4,nanmean(theta(groupIdx{iGroup},1:4),1),'Color',groupColor{iGroup},'LineWidth',3)
%     pause
%
% end


