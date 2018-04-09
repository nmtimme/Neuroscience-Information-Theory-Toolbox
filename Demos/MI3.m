%% Mutual Information Demo 3
% This demo goes through calculating the mutual information for two
% continuous data streams. This example demonstrates the use of predefined
% null models.

% Data Type: Continuous
% Number of Variables: 3
% Information Analysis: Mutual Information
% Trial Based or Single Trial: Single Trial
% Significance Testing: On
% Delays: Off
% Notes: This demo is similar to Demo 7, but a null model from an earlier
%     calculation is used to calculate the p-value.

% We will generate data for three variables. Variables 1 and 2 will be
% related, but the third variable will be independent.
nVar = 3;
nData = 1000;
data = randn([nVar,nData]);
data(2,:) = 0.1*data(2,:) + data(1,:); % variable two is closely related to variable 1, with some noise

% We don't need to use the formattool because the data are already in the
% right format for single trial analysis
DataRaster = data;

% State the data using 4 uniform count bins
MethodAssign = {1,1,'UniCB',{4};...
    1,2,'UniCB',{4};...
    1,3,'UniCB',{4}};
StatesRaster = data2states(DataRaster, MethodAssign);

% Generate the null model using the Monte Carlo outputs from instinfo.
% Because uniform counts binning was used, all MI calculations will use the
% same null model. 
nNullTrials = 10^5;
Method = 'PairMI';
VariableIDs = {1,1,1;1,3,1}; % Variables 1 and 3 (no real interaction)
[waste1,waste2,s,MCInfoVals] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', 'on', 'MCnSamples', nNullTrials, 'MCpThresh', 1);
nullModel = struct;
nullModel.Method = Method;
nullModel.rFact = 10*eps;
nullModel.res = 1/nNullTrials;
nullModel.nCounts = {(nData/4)*ones([4,1]);(nData/4)*ones([4,1])};
MCInfoVals = round(MCInfoVals/nullModel.rFact)*nullModel.rFact; % Round the results to correct for system resolution errors
nullModel.info = unique(MCInfoVals);
nullModel.pmf = histc(MCInfoVals,nullModel.info);

% Perform the information calculation for the single trial across all time
% bins. Use the null model (instinfo tests for applicability of null the
% null model to these data).
Method = 'PairMI';
InfoResults = NaN([2,1]);
p = NaN([2,1]);
s = NaN([2,1]);
VariableIDs = {1,1,1;1,2,1}; % Variables 1 and 2 (real interaction)
[InfoResults(1),p(1),s(1)] = instinfo(StatesRaster, Method, VariableIDs, 'nullModel', nullModel);
VariableIDs = {1,1,1;1,3,1}; % Variables 1 and 3 (no real interaction)
[InfoResults(2),p(2),s(2)] = instinfo(StatesRaster, Method, VariableIDs, 'nullModel', nullModel);

disp(['Variables 1 and 2 (Real Interaction) Mutual Information = ',num2str(InfoResults(1)),' bits (p = ',num2str(p(1)),')'])
disp(['Variables 1 and 3 (No Real Interaction) Mutual Information = ',num2str(InfoResults(2)),' bits (p = ',num2str(p(2)),')'])
