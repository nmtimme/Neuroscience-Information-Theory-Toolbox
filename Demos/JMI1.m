%% Joint Mutual Information Demo 1
% This demo goes through calculating the joint mutual information between
% multiple discrete time series.

% Data Type: Discrete
% Number of Variables: 4
% Information Analysis: Joint Mutual Information
% Trial Based or Single Trial: Single Trial
% Significance Testing: Off
% Delays: On

% For this demo, we will use four variables. We will calculate the mutual
% information between the first two variables jointly and the last two
% variables jointly. We will build a delay into the interaction between the
% variables.
nData = 1000;
data = randi(10,[4,nData]);
copyLocs1 = (rand([1,nData - 1]) < 0.3);
copyLocs2 = (rand([1,nData - 2]) < 0.3);
data(3,[false,copyLocs1]) = data(1,[copyLocs1,false]); % Include a 1 bin delay in interaction between variables 1 and 3
data(4,[false,false,copyLocs2]) = data(2,[copyLocs2,false,false]); % Include a 2 bin delay in interaction between variables 2 and 4

% We don't need to use the formattool because the data are already in the
% right format for single trial analysis
DataRaster = data;

% State the data using 3 uniform width bins
MethodAssign = {1,1,'UniWB',{3};...
    1,2,'UniWB',{3};...
    1,3,'UniWB',{3};...
    1,4,'UniWB',{3}};
StatesRaster = data2states(DataRaster, MethodAssign);

% Perform the information calculation for the single trial across all time
% bins
Method = 'JointMI2';
% Variable 1 and 2 have no delay, but variable 3 and 4 are delayed by 1 and
% 2 time bins. Thus, they have time bin identifiers that are 1 and 2 time
% bins larger than variables 1 and 2. The precise values are not relevant,
% only the relative values to each other.
VariableIDs = {1,1,1;... 
    1,2,1;...
    1,3,2;...
    1,4,3}; 
[InfoResults,p] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', 'on');

disp(['Variables {1,2} and {3,4} (Correct Delays) Joint Mutual Information = ',num2str(InfoResults),' bits (p = ',num2str(p(1)),')'])

% If we rerun the analysis with the wrong delays, the joint mutual
% information collapses
VariableIDs = {1,1,1;... 
    1,2,1;...
    1,3,1;...
    1,4,1}; 
[InfoResults,p] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', 'on');

disp(['Variables {1,2} and {3,4} (Wrong Delays) Joint Mutual Information = ',num2str(InfoResults),' bits (p = ',num2str(p(1)),')'])