%% QUICKTE - quick transfer entropy calculation
% Calculates the transfer entropy between two variables. The data can be
% continuous or discrete. If the data is single trial, the transfer
% entropy is calculated through time. If the data is trial based, the
% transfer entropy is calculated at each time point in the trials. For
% discrete and continuous data, the data will be discretized using uniform
% width bins. The relative time of the past states of the receiver and 
% transmitter can be changed with delays.
%
% Syntax: [transferEntropy,p] = quickTE(dataRec, dataTrans, varargin)
%
% Input:
%   dataRec (double array) - raw data to be analyzed from the receiving
%     variable. The data should be in a two dimensional matrix with time 
%     bins along the first dimension and trials along the second dimension.
%     In other words, each column should contain one trial in time order. 
%     Each row should correspond to an identical time point across trials 
%     (e.g. 5 ms prior to stimulus). The data can be continuous or 
%     discrete. The time and trail data points should match between dataRec
%     and dataTrans. In other words, dataRec(i,j) and dataTrans(i,j) should
%     correspond to measurements taken at the same time point (ith bin) in 
%     the same trial (jth trial).
%   dataTrans (double array) - raw data to be analyzed from the 
%     transmitting variable. The data should be in a two dimensional matrix
%     with time bins along the first dimension and trials along the second 
%     dimension. In other words, each column should contain one trial in 
%     time order. Each row should correspond to an identical time point 
%     across trials (e.g. 5 ms prior to stimulus). The data can be 
%     continuous or discrete. The time and trail data points should match 
%     between dataRec and dataTrans. In other words, dataRec(i,j) and 
%     dataTrans(i,j) should correspond to measurements taken at the same 
%     time point (ith bin) in the same trial (jth trial).
%
% Variable Inputs:
%   (..., 'nbins', nbins) - sets the number uniform width bins used in the
%     discretization (integer, default = 4). 
%   (..., 'delayRec', delayRec) - sets the delay between the past state of
%     the receiver and it's future state (integer, default = 1 (past state
%     immediately precedes the future state)). When delay = 3, the function
%     will compare the future state at time t and the past state at time t
%     - 3. Note, delayRec must be greater than or equal to 1.
%   (..., 'delayTrans', delayTrans) - sets the delay between the past state
%     of the transmitter and the future state of the receiver (integer, 
%     default = 1 (past state immediately preceeds the future state)). When
%     delay = 3, the function will compare the future state of the receiver
%     at time t and the past state of the transmitter at time t - 3. Note, 
%     delayTrans must be greater than or equal to 1.
%   (..., 'MCpThresh', MCpThresh) - sets the p-value cutoff. If the
%     algorithm finds that the p-value will be above the threshold, the
%     calculation ceases (scalar double) (default: 0.001)
%   (..., 'MCnSamples', MCnSamples) - sets the number of Monte Carlo trials
%     to run (scalar double, default: 5000). Sets the resolution of the
%     p-value calculation.
%
% Outputs:
%   transferEntropy (double or double array) - the transfer entropy between 
%     the variables in bits. If trial based data are input, the transfer 
%     entropy will be a vector of length equal to the size of the first
%     dimension of dataRec and dataTrans, where the time bin corresponds to
%     the future state of the receiver. transferEntropy(i) corresponds to 
%     dataY(i,:). Therefore, the first time bin (and possible more 
%     depending on delay settings) will be NaN because the past states of 
%     the receiver and transmitter were unknown. If single trial data are 
%     input, transferEntropy will be a single value representing the 
%     transfer entropy calculated through time.
%   p (double or double array) - the p-value for the corresponding transfer
%     entropy result as measured using Monte Carlo via randomized state
%     observations. A low p-value indicates that the information result is 
%     not likely the result of chance. The randomization proceedure 
%     involved preserves the number of observations, all marginal 
%     probabilities, and the joint probabilities for the receiver states 
%     (the self interaction). The p-value calculation has a resolution of 
%     0.0002 (see MCnSamples above). 
%   
% Examples:
%   Single Trial Normally Distributed Data
%     dataTrans = randn([100,1]);
%     dataRec = [0;dataTrans(1:99) + 0.2*randn([99,1])];
%     [transferEntropy,p] = quickTE(dataRec,dataTrans);
%   Single Trial Normally Distributed Data with a Self Interaction
%     dataRec = repmat(randn([5,1]),[20,1]); % Receiver repeats a pattern
%     dataTrans = [dataRec(2:100) + 0.2*randn([99,1]);0]; % Impose short interaction
%     [transferEntropy,p] = quickTE(dataRec,dataTrans); % delay = 1 produces significant TE
%     [transferEntropy,p] = quickTE(dataRec,dataTrans,'delayRec',5); % delay = 5 reveals TE = 0 by self interaction
%   Trial Based (20 Trials) Uniform Distributed Independent Data
%     dataRec = rand([100,20]);
%     dataTrans = rand([100,20]);
%     [transferEntropy,p] = quickTE(dataRec,dataTrans);
%
% Other m-files required: data2states, instinfo
% Subfunctions: none
% MAT-files required: none
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% May 2016; Last revision: 18-May-2016


function [transferEntropy,p] = quickTE(dataRec,dataTrans,varargin)
%% Parse command line for parameters
nBins = 4;
delayRec = 1;
delayTrans = 1;
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
        case 'delayRec',        delayRec = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'delayTrans',      delayTrans = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'MCpThresh',       MCpThresh = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'MCnSamples',      MCnSamples = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(QUICKTE) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Calculate the Entropy

% Find the number of time bins and trials
nTime = size(dataRec,1);
nTrials = size(dataRec,2);

% Error check the number of time bins and trials
if nTime ~= size(dataTrans,1)
    error('dataRec and dataTrans must have the same number of time bins.')
end
if nTrials ~= size(dataTrans,2)
    error('dataRec and dataTrans must have the same number of trials.')
end
if delayRec < 1
    error('delayRec must be greater than or equal to 1.')
end
if delayTrans < 1
    error('delayTrans must be greater than or equal to 1.')
end
if nTime == 1
    error('The recording has only one time bin, so Transfer Entropy cannot be calculated.')
end
if delayRec >= nTime
    error('delayRec must be less than the number of time bins.')
end
if delayTrans >= nTime
    error('delayTrans must be less than the number of time bins.')
end

% Make a DataRaster
DataRaster = zeros([2,nTime,nTrials]);

% Put the data in the data raster
DataRaster(1,:,:) = dataRec;
DataRaster(2,:,:) = dataTrans;

% State the data using uniform width bins
MethodAssign = {1,1,'UniWB',{nBins};...
    1,2,'UniWB',{nBins}};
StatesRaster = data2states(DataRaster, MethodAssign);
% !!! Note !!! Many other stating options exist. See help data2states.

% Calculate the entropy
Method = 'TE';
if nTrials == 1
    
    % Single Trial
    if delayRec >= delayTrans
        VariableIDs = {1,1,1 + max([delayRec,delayTrans]);...
            1,1,1;...
            1,2,1 + (delayRec - delayTrans)};
    else
        VariableIDs = {1,1,1 + max([delayRec,delayTrans]);...
            1,1,1 + (delayTrans - delayRec);...
            1,2,1};
    end
    [transferEntropy,p] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', MCOpt, 'MCpThresh', MCpThresh, 'MCnSamples', MCnSamples);
    
else
    
    % Trial Based
    transferEntropy = NaN([nTime,1]);
    p = NaN([nTime,1]);
    for iT = (1 + max([delayRec,delayTrans])):nTime
        VariableIDs = {1,1,iT;...
            1,1,iT - delayRec;...
            1,2,iT - delayTrans};
        [transferEntropy(iT),p(iT)] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', MCOpt, 'MCpThresh', MCpThresh, 'MCnSamples', MCnSamples);
    end
    
end










