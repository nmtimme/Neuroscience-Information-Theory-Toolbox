%% QUICKPID - quick partial information decomposition calculation
% Calculates the partial information terms (synergy, redundancy, and unique
% information) between three variables (a Y variable and two X variables).
% Depending on the situation, the information results could be interpretted
% as converging (X variables send information to the Y variable which
% processes the inputs) or as diverging (Y variable sends information that
% is represented in the X variables). The data can be continuous or 
% discrete. If the data is single trial, the partial information terms are 
% calculated through time. If the data is trial based, the partial 
% information terms at each time point in the trials. For discrete and 
% continuous data, the data will be discretized using uniform width bins. 
% The relative time between the states of the variables can be changed with
% delays.
%
% Syntax: [PILVals,p] = quickPID(dataY,dataX1,dataX2,varargin)
%
% Input:
%   dataY (double array) - raw data to be analyzed from the Y variable. The
%     data should be in a two dimensional matrxi with time bins along the
%     first dimension and trials along the second dimension. In other
%     words, each column should contain one trial in time order. Each row
%     should correspond to an identical time point across trials (e.g. 5 ms
%     prior to stimulus). The data can be continuous or discrete. The time
%     and trial data points should match between dataY, dataX1, and dataX2.
%     In other words, dataY(i,j), dataX1(i,j), and dataX2(i,j) should
%     correspond to measurements taken at the same time point (ith bin) in
%     the same trial (jth trial). 
%   dataX1 (double array) - raw data to be analyzed from the X1 variable. 
%     The data should be in a two dimensional matrxi with time bins along 
%     the first dimension and trials along the second dimension. In other
%     words, each column should contain one trial in time order. Each row
%     should correspond to an identical time point across trials (e.g. 5 ms
%     prior to stimulus). The data can be continuous or discrete. The time
%     and trial data points should match between dataY, dataX1, and dataX2.
%     In other words, dataY(i,j), dataX1(i,j), and dataX2(i,j) should
%     correspond to measurements taken at the same time point (ith bin) in
%     the same trial (jth trial). 
%   dataX2 (double array) - raw data to be analyzed from the X2 variable. 
%     The data should be in a two dimensional matrxi with time bins along 
%     the first dimension and trials along the second dimension. In other
%     words, each column should contain one trial in time order. Each row
%     should correspond to an identical time point across trials (e.g. 5 ms
%     prior to stimulus). The data can be continuous or discrete. The time
%     and trial data points should match between dataY, dataX1, and dataX2.
%     In other words, dataY(i,j), dataX1(i,j), and dataX2(i,j) should
%     correspond to measurements taken at the same time point (ith bin) in
%     the same trial (jth trial). 
%
% Variable Inputs:
%   (..., 'nbins', nbins) - sets the number uniform width bins used in the
%     discretization (integer, default = 4). 
%   (..., 'delayX1', delayX1) - sets the delay between the Y state and the 
%     X1 state relative to the Y state (integer, default = 0 (X1 and Y 
%     occur cotemporaneously)). When delay = 3, the function will compare 
%     the Y state at time t and the X1 state at time t - 3. 
%   (..., 'delayX2', delayX2) - sets the delay between the Y state and the 
%     X2 state relative to the Y state (integer, default = 0 (X2 and Y 
%     occur cotemporaneously)). When delay = 3, the function will compare 
%     the Y state at time t and the X2 state at time t - 3. 
%   (..., 'MCpThresh', MCpThresh) - sets the p-value cutoff. If the
%     algorithm finds that the p-value will be above the threshold, the
%     calculation ceases (scalar double) (default: 0.001)
%   (..., 'MCnSamples', MCnSamples) - sets the number of Monte Carlo trials
%     to run (scalar double, default: 5000). Sets the resolution of the
%     p-value calculation.
%
% Outputs:
%   PILVals (double array) - the partial information terms between the 
%     variables in bits. If trial based data are input, the partial 
%     information terms will be an array with first dimension equal to the 
%     size of the first dimension of dataY, dataX1, and dataX2, where the 
%     time bin corresponds to the Y state. Some time bins may have NaN 
%     values depending on delay settings. The second dimension of the array
%     will have 4 elements that correspond to redundancy, unique X1, unique
%     X2, and synergy in that order. If single trial data are input, 
%     PILVals will have four values corresponding to the partial 
%     information terms calculated through time.
%   p (double or double array) - the p-value for the corresponding partial 
%     information result as measured using Monte Carlo via randomized state
%     observations. A low p-value indicates that the information result is 
%     not likely the result of chance. The randomization proceedure 
%     involved preserves the number of observations and all marginal 
%     probabilities. The p-value calculation has a resolution of 
%     0.0002 (see MCnSamples above). 
%   
% Examples:
%   Single Trial Normally Distributed Data (Only X1 Influences Y)
%     dataX1 = randn([100,1]);
%     dataX2 = randn([100,1]);
%     dataY = dataX1 + 0.2*randn(size(dataX1));
%     [PILVals,p] = quickPID(dataY,dataX1,dataX2);
%   Single Trial Redundant Operation
%     dataX1 = randn([100,1]);
%     dataX2 = dataX1;
%     dataY = dataX1 + 0.2*randn(size(dataX1));
%     [PILVals,p] = quickPID(dataY,dataX1,dataX2);
%   Single Trial Synergistic Operation
%     dataX1 = randi(2,[100,1]);
%     dataX2 = randi(2,[100,1]);
%     dataY = dataX1 + dataX2;
%     dataY(dataY == 3) = 5;
%     dataY(dataY < 5) = 1;
%     [PILVals,p] = quickPID(dataY,dataX1,dataX2);
%   Trial Based (20 Trials) Uniform Distributed Independent Data
%     dataY = rand([100,20]);
%     dataX1 = rand([100,20]);
%     dataX2 = rand([100,20]);
%     [PILVals,p] = quickPID(dataY,dataX1,dataX2);
%
% Other m-files required: data2states, instinfo
% Subfunctions: none
% MAT-files required: none
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% May 2016; Last revision: 19-May-2016


function [PILVals,p] = quickPID(dataY,dataX1,dataX2,varargin)
%% Parse command line for parameters
nBins = 4;
delayX1 = 0;
delayX2 = 0;
MCpThresh = 0.001;
MCnSamples = 5000;
if nargout == 2
    MCOpt = 'on';
else
    MCOpt = 'off';
end

iVarArg = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'nBins',           nBins = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'delayX1',         delayX1 = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'delayX2',         delayX2 = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'MCpThresh',       MCpThresh = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'MCnSamples',      MCnSamples = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(QUICKPID) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Calculate the Entropy

% Find the number of time bins and trials
nTime = size(dataY,1);
nTrials = size(dataY,2);

% Error check the number of time bins and trials
if nTime ~= size(dataX1,1)
    error('dataY and dataX1 must have the same number of time bins.')
elseif nTime ~= size(dataX2,1)
    error('dataY and dataX2 must have the same number of time bins.')
end
if nTrials ~= size(dataX1,2)
    error('dataY and dataX1 must have the same number of trials.')
elseif nTrials ~= size(dataX2,2)
    error('dataY and dataX2 must have the same number of trials.')
end
if delayX1 >= nTime
    error('delayX1 must be less than the number of time bins.')
end
if delayX2 >= nTime
    error('delayX2 must be less than the number of time bins.')
end

% Make a DataRaster
DataRaster = zeros([3,nTime,nTrials]);

% Put the data in the data raster
DataRaster(1,:,:) = dataY;
DataRaster(2,:,:) = dataX1;
DataRaster(3,:,:) = dataX2;

% State the data using uniform width bins
MethodAssign = {1,1,'UniWB',{nBins};...
    1,2,'UniWB',{nBins};...
    1,3,'UniWB',{nBins}};
StatesRaster = data2states(DataRaster, MethodAssign);
% !!! Note !!! Many other stating options exist. See help data2states.

% Calculate the entropy
Method = '2PID';
if nTrials == 1
    
    % Single Trial
    if max([delayX1,delayX2]) > 0
        if delayX1 >= delayX2
            VariableIDs = {1,1,1 + max([delayX1,delayX2]);...
                1,2,1;...
                1,3,1 + (delayX1 - delayX2)};
        else
            VariableIDs = {1,1,1 + max([delayX1,delayX2]);...
                1,2,1 + (delayX2 - delayX1);...
                1,3,1};
        end
    else
        VariableIDs = {1,1,1;...
            1,2,1 - delayX1;...
            1,3,1 - delayX2};
    end
    [PILVals,p] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', MCOpt, 'MCpThresh', MCpThresh, 'MCnSamples', MCnSamples);
    PILVals = PILVals';
    p = p';
    
else
    
    % Trial Based
    PILVals = NaN([nTime,4]);
    p = NaN([nTime,4]);
    if max([delayX1,delayX2]) > 0
        tStartOffSet = max([delayX1,delayX2]);
    else
        tStartOffSet = 0;
    end
    if min([delayX1,delayX2]) < 0
        tEndOffSet = min([delayX1,delayX2]);
    else
        tEndOffSet = 0;
    end
    for iT = (1 + tStartOffSet):(nTime + tEndOffSet)
        VariableIDs = {1,1,iT;...
            1,2,iT - delayX1;...
            1,3,iT - delayX2};
        [PILVals(iT,:),p(iT,:)] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', MCOpt, 'MCpThresh', MCpThresh, 'MCnSamples', MCnSamples);
    end
    
end










