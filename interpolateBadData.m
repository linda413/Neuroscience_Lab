function [interpolatedData,badIdx] = interpolateBadData(LFP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [interpolatedData,badIdx] = interpolateBadData(LFP)
%
% Function to handle bad data by interpolation -> it better handles edge 
% effects from fft
% Input
%   LFP - structure containing EEG data
%
% Output
%   interpolatedData - the EEG signal with interpolated sections based on
%       LFP.bad_intervals
%   badIdx - all indices of the EEG singnal that were interpolated
%
% Amber Schedlbauer - 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

badIdx = [];
interpolatedData = LFP.values;

for iBad = 1:size(LFP.bad_intervals,1)
    
    % This removes the bad interval [1,2] that Matt inserts for some reason.
    if LFP.bad_intervals(iBad,1) == 1 && LFP.bad_intervals(iBad,2) == 2
        continue
    end
        
    if LFP.bad_intervals(iBad,1) == 1 || LFP.bad_intervals(iBad,2) == length(LFP.values)
        interpolatedData(LFP.bad_intervals(iBad,1):LFP.bad_intervals(iBad,2)) = 0;
    else
        x = 1:length(LFP.values);
        y = interpolatedData;
        
        x(LFP.bad_intervals(iBad,1):LFP.bad_intervals(iBad,2)) = [];
        y(LFP.bad_intervals(iBad,1):LFP.bad_intervals(iBad,2)) = [];
        
        z = interp1(x,y,1:length(LFP.values),'pchip');
        interpolatedData = z;
    end
    
    badIdx = [badIdx LFP.bad_intervals(iBad,1):LFP.bad_intervals(iBad,2)];
    
end

end