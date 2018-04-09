%% Entropy Demo 1
% This demo goes through calculating the entropy for time locked trials
% from a single continuous time series. 

% Data Type: Continuous
% Number of Variables: 1
% Information Analysis: Entropy
% Trial Based or Single Trial: Trial Based
% Significance Testing: N/A
% Delays: N/A

% Generate the data
nData = 10000;
data = randn([1,nData]);
time = 1:nData; % assume the data are binned at 1 time unit (could be milliseconds)
timelocks = 100:100:(nData - 100); % assume each trial was conducted every 100 time units
binsize = 1; % leave the bins unchanged
maxLead = 20; % assume the trial started 20 time units before the time lock point
maxLag = 30; % assume the trial ended 30 time units after the time lock point
method = 'ave'; % the method is not relevant because we did not change the bin size

% Format the data based on the trials
[DataRaster,timeboundaries] = formattool(data, time, timelocks, binsize, maxLead, maxLag, method);

% State the data using 4 uniform width bins
MethodAssign = {1,1,'UniWB',{4}};
StatesRaster = data2states(DataRaster, MethodAssign);

% Perform the information calculation for all time bins
Method = 'Ent';
nT = size(StatesRaster,2);
InfoResults = NaN([nT,1]);
for iT = 1:nT
    VariableIDs = {1,1,iT};
    InfoResults(iT) = instinfo(StatesRaster, Method, VariableIDs);
end

figure
plot(mean(timeboundaries,2),InfoResults)
xlabel('Time')
ylabel('Information (bits)')
title('Entropy')