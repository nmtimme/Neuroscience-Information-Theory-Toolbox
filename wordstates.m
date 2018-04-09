%% WORDSTATES - convert single time bin states to words
% Converts single time bin states to words (joint vector valued state
% combinations). For instance, the state time sequence [1,1,2] is
% distinguished from [1,2,2]. Note, wordstates will reduce the total number
% of time bins in the StatesRaster. Trial based data will be converted
% across trials at each individual time. Single trial data will be
% converted through time.
%
% Syntax: [StatesRaster, timeboundaries] = wordstates(StatesRaster, WordInfo, varargin)
%
% Input:
%   StatesRaster (cell array or double array) - trial state data. If a
%     double array is used, it should be number of variables by number of
%     time bins by number of trials. Each element should be an integer
%     state number (state number = 1, 2, 3, ...). If a cell array is used,
%     it should have only one dimension and each element should be a double
%     array with the dimensions listed above. Each element of the cell
%     array is referred to as a 'data category'.
%   WordInfo (integer array) - sets the number of single time bin states to
%     combine into words. The array should have dimensions number of data
%     categories by 2. Each row controls the word conversion for a data
%     category. The first column is the data category and the second
%     category is the number of contiguous single time bin states to
%     combine. If StatesRaster contains only one data category, the first
%     column of WordInfo can be removed. Note, all variables in a given
%     data category will be converted. Also, note that data categories not
%     listed in WordInfo will not undergo word conversion.
%
% Variable Inputs:
%   (..., timeboundaries) - records the time boundaries for the time bins.
%     This must be a cell array or double array (matching StatesRaster).
%     Each element (in the case of a cell array) should have dimension
%     number of time bins (matching the data category from StatesRaster) by
%     2. 
%
% Outputs:
%   StatesRaster (cell array or double array) - trial state data after word
%     conversion. It will be the same data type as the StatesRaster input,
%     though the time dimension (second dimension) of the double array will
%     be smaller for the data categories that underwent word conversion.
%   timeboundaries (cell array or double array) - new time boundaries
%     following word conversion. If timeboundaries was not input, this
%     variable will be NaN. It will have the same structure as the input
%     timeboundaries, but it will have fewer bins for the data categories
%     that underwent word conversion.
%   
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% May 2016; Last revision: 13-May-2016


function [StatesRaster, timeboundaries] = wordstates(StatesRaster, WordInfo, varargin)
%% Parse command line for parameters
TBOpt = false;
timeboundaries = NaN;

if length(varargin) == 1
    timeboundaries = varargin{1};
    TBOpt = true;
end


%% Perform Initial Operations

% Error check inputs and possible reorganize to ease later processing
if iscell(StatesRaster)
    Back2ArrayFlag = false;
    if length(StatesRaster) ~= 1
        if size(WordInfo,2) ~= 2
            error('WordInfo is not the appropriate size. Probably the data category is missing.')
        end
    else
        if size(WordInfo,2) == 1
            WordInfo = [ones(size(WordInfo)),WordInfo];
        elseif size(WordInfo,2) ~= 2
            error('WordInfo is not the appropriate size.')
        end
    end
    if TBOpt
        if length(timeboundaries) ~= length(StatesRaster)
            error('timeboundaries and StatesRasters do not have the same number of data categories.')
        else
            for iDC = 1:length(StatesRaster)
                if size(timeboundaries{iDC},1) ~= size(StatesRaster{iDC},2)
                    error(['timeboundaries and StatesRaster are inconsistent in data category ',num2str(iDC)])
                end
            end
        end
    end
else
    Back2ArrayFlag = true;
    StatesRaster = {StatesRaster};
    if size(WordInfo,2) == 1
        WordInfo = [ones(size(WordInfo)),WordInfo];
    end
    if TBOpt
        if size(timeboundaries,1) ~= size(StatesRaster{1},2)
            error('timeboundaries and StatesRaster are inconsistent')
        end
    end
end

% Make sure we didn't get duplicate instructions
if length(unique(WordInfo(:,1))) ~= length(WordInfo(:,1))
    error('Duplicate data category instructions in WordInfo.')
end



%% Perform the Word Conversion

for iConv = 1:size(WordInfo,1)
    
    if size(StatesRaster{WordInfo(iConv,1)},3) > 1 % Trial Based Data
        
        % Combine states for unique words across trials
        nT = size(StatesRaster{WordInfo(iConv,1)},2);
        Temp1 = NaN(size(StatesRaster{WordInfo(iConv,1)}) - [0,WordInfo(iConv,2) - 1,0]);
        for iVar = 1:size(StatesRaster{WordInfo(iConv,1)},1)
            for iT = 1:(nT - WordInfo(iConv,2) + 1)
                Temp2 = squeeze(StatesRaster{WordInfo(iConv,1)}(iVar,iT:(iT + WordInfo(iConv,2) - 1),:))';
                [B,I,Temp1(iVar,iT,:)] = unique(Temp2,'rows');
            end
        end
        StatesRaster{WordInfo(iConv,1)} = Temp1;
        
    else % Single Trial Data
        
        % Combine states for unique words across time
        nT = size(StatesRaster{WordInfo(iConv,1)},2);
        Temp1 = NaN(size(StatesRaster{WordInfo(iConv,1)}) - [0,WordInfo(iConv,2) - 1]);
        for iVar = 1:size(StatesRaster{WordInfo(iConv,1)},1)
            Temp2 = NaN([nT - WordInfo(iConv,2) + 1,WordInfo(iConv,2)]);
            for iT = 1:WordInfo(iConv,2)
                Temp2(:,iT) = StatesRaster{WordInfo(iConv,1)}(iVar,iT:(nT - WordInfo(iConv,2) + iT));
            end
            [B,I,Temp1(iVar,:)] = unique(Temp2,'rows');
        end
        StatesRaster{WordInfo(iConv,1)} = Temp1;
        
    end
    
    if TBOpt % Correct the time boundaries
        
        timeboundaries{WordInfo(iConv,1)}(1:(nT - WordInfo(iConv,2) + 1),2) = timeboundaries{WordInfo(iConv,1)}(WordInfo(iConv,2):nT,2);
        timeboundaries{WordInfo(iConv,1)}((nT - WordInfo(iConv,2) + 2):nT,:) = [];
        
    end
    
end

% Convert StatesRaster back to an array, if necessary
if Back2ArrayFlag
    StatesRaster = StatesRaster{1};
end

