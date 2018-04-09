%% Mutual Information Demo 2
% This demo goes through calculating the mutual information for three 
% continuous time series.

% Data Type: Continuous
% Number of Variables: 3
% Information Analysis: Mutual Information
% Trial Based or Single Trial: Single Trial
% Significance Testing: On
% Delays: Off

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

% Perform the information calculation for the single trial across all time
% bins
Method = 'PairMI';
InfoResults = NaN([2,1]);
p = NaN([2,1]);
VariableIDs = {1,1,1;1,2,1}; % Variables 1 and 2 (real interaction)
[InfoResults(1),p(1)] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', 'on');
VariableIDs = {1,1,1;1,3,1}; % Variables 1 and 3 (no real interaction)
[InfoResults(2),p(2)] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', 'on');

disp(['Variables 1 and 2 (Real Interaction) Mutual Information = ',num2str(InfoResults(1)),' bits (p = ',num2str(p(1)),')'])
disp(['Variables 1 and 3 (No Real Interaction) Mutual Information = ',num2str(InfoResults(2)),' bits (p = ',num2str(p(2)),')'])