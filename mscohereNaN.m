function [Cxy,F] = mscohereNaN(x,y,hamWin,overlap,nFFT,sFreq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Cxy,F] = mscohereNaN(x,y,hamWin,overlap,nFFT,sFreq)
%
% Inputs like mscohere, but the nans will be ignored unless they are more
% than 50% of the usable recording. Segments too short are defined as being
% twice the Hamming window in length. They are ignored. Whenever an input 
% is unusable, mscohereNaN outputs a nan.
%
% Matthew B Schmit (Mattenator) 2016 from Cowen's code
% Edited by Amber Schedlbauer 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for nans
IX = isnan(x+y);
if mean(IX) > 0.5
    warning('Window more that 50% NaN!')
    Cxy = nan(nFFT/2+1,1);
end

% Define minimum length for calculation
minLength = hamWin*2;

if ~any(IX) % no nans
    [Cxy,F] = mscohere(x,y,hamWin,overlap,nFFT,sFreq);  
else
    
    % Finds nonNaN intervals so the coherence can be calculated
    intervals = cat(2,find((-diff([1 IX]) > 0.5))',find(diff([IX 1]) > 0.5)');
    
    % Coherence for that data segment that does not contain NaNs
    ctr = 1;
    for iSeg = 1:size(intervals,1)
        
        % Size of analyzeable data in interval
        subWindowSize(iSeg) = intervals(iSeg,2) - intervals(iSeg,1);
        
        if subWindowSize(iSeg) > minLength
            [Cxy(:,ctr),F] = mscohere(x(intervals(iSeg,1):intervals(iSeg,2)),y(intervals(iSeg,1):intervals(iSeg,2)),hamWin,overlap,nFFT,sFreq);
            ctr = ctr + 1;
        end
        
    end
   
    subWindowSize(subWindowSize < minLength) = [];
    
    if sum(subWindowSize)/length(x) < 0.5
        display('Window does not have enough long contiguous chunks.')
        Cxy = nan(nFFT/2+1,1);
    else
        subWindowSize = subWindowSize/sum(subWindowSize);
        Cxy = sum(Cxy.*repmat(subWindowSize,size(Cxy,1),1),2);
    end
    
end

end
    
    