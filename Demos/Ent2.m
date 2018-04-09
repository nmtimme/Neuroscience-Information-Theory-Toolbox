%% Entropy Demo 2
% This demo goes through calculating the entropy using a single discrete 
% time series.

% Data Type: Discrete
% Number of Variables: 1
% Information Analysis: Entropy
% Trial Based or Single Trial: Single Trial
% Significance Testing: N/A
% Delays: N/A

% Generate the data
nData = 100;
data = randi(10,[1,nData]);
time = 1:nData; % assume the data are binned at 1 time unit (could be milliseconds)

% We don't need to use the formattool because the data are already in the
% right format for single trial analysis
DataRaster = data;

% State the data using 4 uniform count bins
MethodAssign = {1,1,'UniCB',{4}};
StatesRaster = data2states(DataRaster, MethodAssign);

% Perform the information calculation for the single trial across all time
% bins
Method = 'Ent';
VariableIDs = {1,1,1};
InfoResults = instinfo(StatesRaster, Method, VariableIDs);

disp(['The entropy of the data was ',num2str(InfoResults),' bits.'])