function UDIs = plotAnalyses(dataCell,frequencies,UDIs,plotType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The important thing is to get the data in the correct format!
% You should really plot no more than 4 groups for the matrices.
%
% UDIs = structure with plotting inputs
%   UDIs.groups = {};
%   UDIs.groupColor = {};
%   UDIs.XAXIS = [];
%   UDIs.YAXIS = [];
%   UDIs.COLORMAP = '';
%   UDIs.XLABEL = {''};
%   UDIs.YLABEL = '';
%   UDIs.TITLE = {''};
%   UDIs.COLORBARTITLE = '';
%   UDIs.SAVENAME = '';
%   Note: Pay attention to the data type for each field!
%
% plotType: 'boxPlotScatter' 'matrix' 'PSD' 'PSDsig' 'histogramPlot'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extra tidbits

% Used to draw dotted lines on the plots to highlight the theta band
frequenciesOfInterest = [6 10];
[~,thetaBand(1)] = min(abs(frequencies - frequenciesOfInterest(1)));
[~,thetaBand(2)] = min(abs(frequencies - frequenciesOfInterest(2)));

%% Plotting

switch plotType
    
    case 'boxPlotScatter' 
        
        %%% Single Box Plot with Scatter Points %%%
        % dataCell = cell array with each cell containing a n x 1 (animals x 1) matrix of the data for that group
        % UDIs needed: groups, groupColor, SAVENAME
        
        figure('Position',[100 100 500 500])
        plotBoxScatter(dataCell,UDIs.groups,UDIs.groupColor)
        set(gca,'Position',[0.13 0.12 0.8 0.8])
        set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
        xlabel(UDIs.XLABEL{1})
        ylabel(UDIs.YLABEL)
        title(UDIs.TITLE{1})
        box off
        save2pdf(UDIs.SAVENAME)
        %%%
        
    case 'matrix'
        
        %%% Matrix (spectrograms, P-episode matrices, coherence matrices) %%%
        % dataCell = cell array with each cell containing n x m (freq x time) matrix of the data for that group
        % UDIs needed: XAXIS, YAXIS, COLORMAP, SAVENAME
        
        figure('Position',[100 100 560*length(dataCell) 500]);
        colormap(UDIs.COLORMAP)
        
        for iGroup = 1:length(dataCell)
            
            handleAX(iGroup) = subplot(1,length(dataCell),iGroup);
            %contourf(UDIs.XAXIS,1:length(UDIs.YAXIS),dataCell{iGroup},30,'LineColor','none')
            imagesc(UDIs.XAXIS,1:length(UDIs.YAXIS),dataCell{iGroup})
            hline(thetaBand,{'--k','--k'})
            axis xy
            set(gca,'YTick',1:5:length(UDIs.YAXIS),'YTickLabel',[])
            caxis([min(cellfun(@(x) min(x(:)),dataCell)) max(cellfun(@(x) max(x(:)),dataCell))])
            
            if iGroup == 1
                set(gca,'YTickLabel',round(UDIs.YAXIS(1:5:end),1))
                ylabel(UDIs.YLABEL)
            end
            
            if iGroup == length(dataCell)
                handleCOLORBAR = colorbar;
                title(handleCOLORBAR,UDIs.COLORBARTITLE)
            end
            
            handleAX(iGroup).Position = [handleAX(iGroup).Position(1) 0.12 (400/(560*length(dataCell))) 0.8];
            xlabel(UDIs.XLABEL{iGroup})
            
            set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
            title(UDIs.TITLE{iGroup})
        end
        
        save2pdf(UDIs.SAVENAME)
        %%%
        
    case 'PSD'
        
        %%% Single plot PSD %%%
        % dataCell = cell array with each cell containing a n x 1 (freq x 1) matrix of the data for that group
        % UDIs needed: groups, groupColor, XAXIS, SAVENAME
        %
        % Ex: PSD = squeeze(nanmean(powerAllConvert,2));
        %     dataCell{1} = nanmean(PSD(:,groupIdx{1}),2); dataCell{2} = nanmean(PSD(:,groupIdx{2}),2);
        
        figure('Position',[100 100 500 500])
        hold on
        for iGroup = 1:length(dataCell)
            plot(UDIs.XAXIS,dataCell{iGroup},'LineWidth',6,'Color',UDIs.groupColor{iGroup})
        end
        vline(frequenciesOfInterest,'k--')
        set(gca,'Position',[0.13 0.12 0.8 0.8])
        set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
        xlabel(UDIs.XLABEL{1})
        ylabel(UDIs.YLABEL)
        legend(UDIs.groups)
        title(UDIs.TITLE{1})
        save2pdf(UDIs.SAVENAME)
        %%%
        
    case 'PSDsig'
        
        %%% Single plot PSD with Significance Shading - 2 groups ONLY! %%%
        % dataCell = cell array with each cell containing a n x m (animal x freq) matrix of the data for that group
        % UDIs needed: groups, groupColor, XAXIS, SAVENAME
        %
        % Ex: PSD = squeeze(nanmean(powerAllConvert,2))';
        %     dataCell{1} = PSD(spikeGroupIdx{1},:); dataCell{2} = PSD(spikeGroupIdx{2},:);
        
        % Uses a Welch's ttest for only the frequencies in the theta range and corrects for multiple comparisons using Bonferroni
        [hWeird,~] = ttest2(dataCell{1},dataCell{2},'Vartype','unequal','Alpha',0.05/numel(thetaBand(1):thetaBand(2)));
        hWeird(isnan(hWeird)) = 0;
        freqWin = zeros(1,length(UDIs.XAXIS)); freqWin(thetaBand(1):thetaBand(2)) = 1; hWeird = logical(hWeird.*freqWin);

        dataCell = cellfun(@(x) nanmean(x,1)',dataCell,'UniformOutput',false);

        figure('Position',[100 100 500 500]);
        hold on
        
        if sum(hWeird) == 1
            plot([UDIs.XAXIS(hWeird) UDIs.XAXIS(hWeird)],[dataCell{1}(hWeird) flip(dataCell{2}(hWeird))],'Color',[0.5 0.5 0.5 0.15],'LineWidth',5);
        else
            fill([UDIs.XAXIS(hWeird) flip(UDIs.XAXIS(hWeird))],[dataCell{1}(hWeird)' flip(dataCell{2}(hWeird))'],'k','LineStyle','none')
        end
        
        alpha(0.15)
        
        handleLINE = plot(UDIs.XAXIS,[dataCell{1} dataCell{2}],'LineWidth',4);
        handleLINE(1).Color = UDIs.groupColor{1}; handleLINE(2).Color = UDIs.groupColor{2};
        vline(frequenciesOfInterest,'k--')
        
        set(gca,'Position',[0.13 0.12 0.8 0.8])
        set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
        xlabel(UDIs.XLABEL{1})
        ylabel(UDIs.YLABEL)
        legend(handleLINE,UDIs.groups)
        title(UDIs.TITLE{1})
        save2pdf(UDIs.SAVENAME)
        %%%
        
        %%% PSD with Subplot of T-statistic
        % Place a "keyboard" on Line 116
%         [~,~,~,stats] = ttest2(dataCell{1},dataCell{2},'Vartype','unequal');
%         
%         figure('Position',[100 100 500 500])
%         hold on
%         
%         plot(frequencies,stats.tstat,'k','LineWidth',4)
%         degreesFreedom = min(round(stats.df));
%         hline([tinv(0.025,degreesFreedom) tinv(0.975,degreesFreedom)])
%         
%         set(gca,'Position',[0.13 0.12 0.8 0.8])
%         set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
%         ylabel('t statistic')
%         xlabel('Frequencies (Hz)')
        
        %%%
        
    case 'histogramPlot'
        
        %%% Single plot Histogram %%%
        % dataCell = cell array
        % UDIs needed: groups, XAXIS, SAVENAME
        
        figure('Position',[100 100 1500 500])
        hold on
        for iGroup = 1:length(dataCell)
            histogram(dataCell{iGroup}(:),'FaceColor',UDIs.groupColor{iGroup})%,'FaceAlpha',0.5)
        end
        set(gca,'XTick',1:length(UDIs.XAXIS),'XTickLabel',round(UDIs.XAXIS,1))
        set(gca,'FontName','Helvetica','FontSize',20,'FontWeight','bold')
        xlabel(UDIs.XLABEL{1})
        ylabel(UDIs.YLABEL)
        legend(UDIs.groups)
        title(UDIs.TITLE{1})
        save2pdf(UDIs.SAVENAME)
        %%%
        
end

%% Reset UDIs structure

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

end