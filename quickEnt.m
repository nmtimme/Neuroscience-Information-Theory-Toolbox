%% QUICKENT - quick entropy calculation
% Calculates the entropy for the input data. The data can be continuous or
% discrete. If the data is single trial, the entropy is calculated through
% time. If the data is trial based, the entropy is calculated at each time
% point in the trials. For discrete and continuous data, the data will be
% discretized using uniform width bins.
%
% Syntax: [entropy] = quickEnt(data, varargin)
%
% Input:
%   data (double array) - raw data to be analyzed. The data should be in a
%     two dimensional matrix with time bins along the first dimension and
%     trials along the second dimension. In other words, each column should
%     contain one trial in time order. Each row should correspond to an
%     identical time point across trials (e.g. 5 ms prior to stimulus). The
%     data can be continuous or discrete.
%
% Variable Inputs:
%   (..., 'nbins', nbins) - sets the number uniform width bins used in the
%     discretization (integer, default = 4). 
%
% Outputs:
%   entropy (double or double array) - the entropy of the data in bits. If
%     trial based data are input, the entropy will be a vector of length
%     equal to the size of the first dimension of data.
%   
% Examples:
%   Single Trial Normally Distributed Data
%     data = randn([100,1]);
%     entropy = quickent(data);
%   Trial Based (20 Trials) Uniform Distributed Data
%     data = rand([100,20]);
%     entropy = quickent(data);
%
% Other m-files required: data2states, instinfo
% Subfunctions: none
% MAT-files required: none
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% May 2016; Last revision: 13-May-2016


function [entropy] = quickEnt(data, varargin)
%% Parse command line for parameters
nBins = 4;

iVarArg = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'nBins',        nBins = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
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
nTime = size(data,1);
nTrials = size(data,2);

% Make a DataRaster
DataRaster = zeros([1,nTime,nTrials]);

% Put the data in the data raster
DataRaster(1,:,:) = data;

% State the data using uniform width bins
MethodAssign = {1,1,'UniWB',{nBins}};
StatesRaster = data2states(DataRaster, MethodAssign);
% !!! Note !!! Many other stating options exist. See help data2states.

% Calculate the entropy
Method = 'Ent';
if nTrials == 1
    
    % Single Trial
    VariableIDs = {1,1,1};
    entropy = instinfo(StatesRaster, Method, VariableIDs);
    
else
    
    % Trial Based
    entropy = NaN([nTime,1]);
    for iT = 1:nTime
        VariableIDs = {1,1,iT};
        entropy(iT) = instinfo(StatesRaster, Method, VariableIDs);
    end
    
end










