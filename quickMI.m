%% QUICKMI - quick mutual information calculation
% Calculates the mutual information between two variables. The data can be
% continuous or discrete. If the data is single trial, the mutual
% information is calculated through time. If the data is trial based, the
% mutual information is calculated at each time point in the trials. For
% discrete and continuous data, the data will be discretized using uniform
% width bins. The relative times of the states from the variables can be
% changed with a delay.
%
% Syntax: [mutualInfo,p] = quickMI(data1, data2, varargin)
%
% Input:
%   data1 (double array) - raw data to be analyzed from the first variable.
%     The data should be in a two dimensional matrix with time bins along 
%     the first dimension and trials along the second dimension. In other 
%     words, each column should contain one trial in time order. Each row 
%     should correspond to an identical time point across trials (e.g. 5 ms
%     prior to stimulus). The data can be continuous or discrete. The time 
%     and trial data points should match between data1 and data2. In other
%     words, data1(i,j) and data2(i,j) should correspond to measurements
%     taken at the same time point (ith bin) in the same trial (jth trial).
%   data2 (double array) - raw data to be analyzed from the second 
%     variable. The data should be in a two dimensional matrix with time 
%     bins along the first dimension and trials along the second dimension.
%     In other words, each column should contain one trial in time order. 
%     Each row should correspond to an identical time point across trials 
%     (e.g. 5 ms prior to stimulus). The data can be continuous or 
%     discrete. The time and trail data points should match between data1 
%     and data2. In other words, data1(i,j) and data2(i,j) should 
%     correspond to measurements taken at the same time point (ith bin) in 
%     the same trial (jth trial).
%
% Variable Inputs:
%   (..., 'nbins', nbins) - sets the number uniform width bins used in the
%     discretization (integer, default = 4). 
%   (..., 'delay', delay) - sets the delay for data2 relative to data1 in
%     units of time bins (integer, default = 0). When delay = 1, the
%     function will compare data1 at time t to data2 at time t - 1. 
%   (..., 'MCpThresh', MCpThresh) - sets the p-value cutoff. If the
%     algorithm finds that the p-value will be above the threshold, the
%     calculation ceases (scalar double) (default: 0.001)
%   (..., 'MCnSamples', MCnSamples) - sets the number of Monte Carlo trials
%     to run (scalar double, default: 5000). Sets the resolution of the
%     p-value calculation.
%
% Outputs:
%   mutualInfo (double or double array) - the mutual information between 
%     the variables in bits. If trial based data are input, the mutual 
%     information will be a vector of length equal to the size of the first
%     dimension of data1 and data2. mutualInfo(i) corresponds to data1(i,:)
%     and data2(i - delay,:). If single trial data is input, mutualInfo
%     will be a single value representing the mutual information calculated
%     throughout time.
%   p (double or double array) - the p-value for the corresponding mutual 
%     information result as measured using Monte Carlo via randomized state
%     observations. A low p-value indicates that the information result is 
%     not likely the result of chance. The randomization proceedure 
%     involved preserves the number of observations and the marginal 
%     probabilities. The p-value calculation has a resolution of 0.0002 
%     (see MCnSamples above).
%   
% Examples:
%   Single Trial Normally Distributed Data
%     data1 = randn([100,1]);
%     data2 = data1 + 0.2*randn(size(data1));
%     [mutualInfo,p] = quickMI(data1,data2);
%   Trial Based (20 Trials) Uniform Distributed Independent Data
%     data1 = rand([100,20]);
%     data2 = rand([100,20]);
%     [mutualInfo,p] = quickMI(data1,data2);
%
% Other m-files required: data2states, instinfo
% Subfunctions: none
% MAT-files required: none
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% May 2016; Last revision: 18-May-2016


function [mutualInfo,p] = quickMI(data1,data2,varargin)
%% Parse command line for parameters
nBins = 4;
delay = 0;
if nargout == 2
    MCOpt = 'on';
else
    MCOpt = 'off';
end

iVarArg = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'nBins',        nBins = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'delay',        delay = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'MCpThresh',       MCpThresh = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'MCnSamples',      MCnSamples = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(QUICKENT) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Calculate the Entropy

% Find the number of time bins and trials
nTime = size(data1,1);
nTrials = size(data1,2);

% Error check the number of time bins and trials
if nTime ~= size(data2,1)
    error('data1 and data2 must have the same number of time bins.')
end
if nTrials ~= size(data2,2)
    error('data1 and data2 must have the same number of trials.')
end
if abs(delay) >= nTime
    error('abs(delay) must be less than the number of time bins')
end

% Make a DataRaster
DataRaster = zeros([2,nTime,nTrials]);

% Put the data in the data raster
DataRaster(1,:,:) = data1;
DataRaster(2,:,:) = data2;

% State the data using uniform width bins
MethodAssign = {1,1,'UniWB',{nBins};...
    1,2,'UniWB',{nBins}};
StatesRaster = data2states(DataRaster, MethodAssign);
% !!! Note !!! Many other stating options exist. See help data2states.

% Calculate the entropy
Method = 'PairMI';
if nTrials == 1
    
    % Single Trial
    if delay >= 0
        VariableIDs = {1,1,1 + delay;...
            1,2,1};
    else
        VariableIDs = {1,1,1;...
            1,2,1 - delay};
    end
    [mutualInfo,p] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', MCOpt);
    
else
    
    % Trial Based
    mutualInfo = NaN([nTime,1]);
    p = NaN([nTime,1]);
    if delay >= 0
        for iT = (1 + delay):nTime
            VariableIDs = {1,1,iT;...
                1,2,iT - delay};
            [mutualInfo(iT),p(iT)] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', MCOpt, 'MCpThresh', MCpThresh, 'MCnSamples', MCnSamples);
        end
    else
        for iT = 1:(nTime + delay)
            VariableIDs = {1,1,iT;...
                1,2,iT - delay};
            [mutualInfo(iT),p(iT)] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', MCOpt, 'MCpThresh', MCpThresh, 'MCnSamples', MCnSamples);
        end
    end
    
end










