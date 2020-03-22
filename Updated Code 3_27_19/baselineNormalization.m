function [powerDB,powerZ] = baselineNormalization(powerEntireSpec,powerInterest,analysisAnimal,task)

%%% Stuff is hard-coded based upon specific analyses!!! %%%
%%% Must modify baselineWindow if use Block time window section!!!

%error('Please open up the baselineNormalization.m script and choose the block of code you want to execute (block time window or time windows of specific behaviors).')

sRate = 1000;

%% Block time window (e.g. baselineWindow = [0 300];)
timeStamps = 0:(1/sRate):(size(powerEntireSpec,2)/sRate);
baselineWindow = [0 300]; % 1st 5 minutes of recording session
timeVec = timeStamps >= baselineWindow(1) & timeStamps < baselineWindow(2); % Baseline index
timeVec = logical(timeVec);

%% Time windows of specific behavior

% % Read in data from NOR behavior excel file
% [~,txt,~] = xlsread('/Users/amberschedlbauer/Dropbox/Gurkoff Lab/Behavior/NOR Video Analysis.xlsx',1,'C3:DS102');
% 
% % Extract animal and task listed in the file
% [animalNames,remain] = strtok(txt(:,1),'_');
% remain = cellfun(@(x) x(2:end),remain,'UniformOutput',false);
% [taskNames,~] = strtok(remain,'_');
% 
% % Find the row of the animals of interest
% animalIdx = find(ismember(animalNames,analysisAnimal));
% 
% % Find the row of the task of interest
% taskIdx = find(ismember(taskNames,task));
% 
% % Find the row of both the animal and task of interest
% row = intersect(animalIdx,taskIdx);
% 
% % Behavior used for background window
% behavior = {'NM';'S'};
% 
% if isempty(row)
%     timeWindow = [];
% elseif length(row) == 1
%     
%     tempIdx = [];
%     
%     for iBehavior = 1:length(behavior)
%         
%         if ismember(behavior{iBehavior},{'M';'S'})
%             temp = find(strcmp(txt(row,2:end),behavior{iBehavior}));
%             tempIdx = [tempIdx temp];
%             clear temp
%         else
%             temp = strfind(txt(row,2:end),behavior{iBehavior});
%             tempIdx = [tempIdx find(~cellfun(@isempty,temp))];
%         end
%         
%     end
%     
%     behaviorIdx= sort(tempIdx);
%     
%     % Finds the start and end times of the windows of the behavior
%     if isempty(behaviorIdx)
%         timeWindow = [];
%         endTime = [];
%     else
%         startTime = behaviorIdx([true (diff(behaviorIdx) ~= 1)])';
%         endTime = behaviorIdx([(diff(behaviorIdx) ~= 1) true])';
%         
%         timeWindow = [startTime endTime];
%         
%         % Checks for beginning and end of video analysis
%         if sum(startTime == 120)
%             timeWindow((startTime == 120),:) = [];
%         end
%         
%         if sum(timeWindow(:,2) == 120)
%             timeWindow((timeWindow(:,2) == 120),2) = 119;
%         end
%     end
%     
% end
% 
% timeVec = zeros(1,length(powerEntireSpec));
% for iTime = 1:size(timeWindow,1)
%     timeVec(timeWindow(iTime,1)*sRate:(timeWindow(iTime,2)*sRate + (sRate-1))) = 1;
% end
% timeVec = logical(timeVec);

%% Baseline Normalizations

if mean(isnan(powerEntireSpec(1,timeVec))) > 0.5
    warning('Baseline window is more than 50% artifact!')
end

baselineMean = nanmean(powerEntireSpec(:,timeVec),2);
baselineSTD = nanstd(powerEntireSpec(:,timeVec),[],2);

% Decibel conversion
powerDB = 10*log10(bsxfun(@rdivide,powerInterest,baselineMean));

% Z-transform
powerZ = bsxfun(@minus,powerInterest,baselineMean);
powerZ = bsxfun(@rdivide,powerZ,baselineSTD);

end

