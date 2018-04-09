%% Joint Entropy Demo 1
% This demo goes through calculating the joint entropy for 3 discrete time 
% series.

% Data Type: Discrete
% Number of Variables: 3
% Information Analysis: Joint Entropy
% Trial Based or Single Trial: Single Trial
% Significance Testing: N/A
% Delays: Off

% Generate the data
nData = 100;
data = randi(10,[3,nData]);
time = 1:nData; % assume the data are binned at 1 time unit (could be milliseconds)

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
Method = 'JointEnt';
VariableIDs = {1,1,1;...
    1,2,1;...
    1,3,1};
InfoResults = instinfo(StatesRaster, Method, VariableIDs);

disp(['The joint entropy of the data was ',num2str(InfoResults),' bits.'])