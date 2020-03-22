function plotBoxScatter(dataCell,groups,scatterColors)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotBoxScatter(dataCell,groups,scatterColors)
%
% Plots a box plot superimposed on a scatter plot that shows the
% distribution of data within each group
%
% INPUTS:
%   dataCell - cell array with each group's data in each cell
%   groups - cell array of group names
%   scatterColors - cell array with each cell containing the color each
%   group should be plotted
%
% Amber Schedlbauer - 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If scatterColors is empty, plots the points in black (default).
if isempty(scatterColors)
    scatterColors = cell(1,length(dataCell));
    scatterColors(:) = {'k'};
end

% If the groups contain different number of observations, nan pads the data
% to be able to use the boxplot function.
maxSize = max(cellfun(@numel,dataCell));    % Get the maximum vector size
fcn = @(x) [x nan(1,maxSize-numel(x))];  % Create an anonymous function
rmat = cellfun(fcn,dataCell,'UniformOutput',false);  % Pad each cell with NaNs
rmat = vertcat(rmat{:})';

% Create figure
%figure('Position',[100 100 500 500])
%hold on

% Plots the individual data points
for iPlot = 1:length(dataCell)
    scatter(ones(1,length(dataCell{iPlot}))*iPlot,dataCell{iPlot},50,'filled','MarkerFaceColor',scatterColors{iPlot})
    hold on
end

% Plots the box plots for the groups
boxplot(rmat,groups)

end