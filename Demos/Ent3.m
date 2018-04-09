%% Entropy Demo 3
% This demo goes through calculating the entropy for time locked trials
% from multiple continuous time series with multiple user defined 
% categories. 

% Data Type: Continuous
% Number of Variables: 6 (4 of one category, 2 of another)
% Information Analysis: Entropy
% Trial Based or Single Trial: Trial Based
% Significance Testing: N/A
% Delays: N/A

% Generate the data
nVar = 6;
nData = 10000;
data = randn([nVar,nData]);
time = 1:nData; % assume the data are binned at 1 time unit (could be milliseconds)
timelocks = 100:100:(nData - 100); % assume each trial was conducted every 100 time units
binsize = 1; % leave the bins unchanged
maxLead = 20; % assume the trial started 20 time units before the time lock point
maxLag = 30; % assume the trial ended 30 time units after the time lock point
method = 'count'; % the method is not relevant because we did not change the bin size

% Format the data based on the trials
DataRaster = cell([2,1]); % there are two categories of data
timeboundaries = cell([2,1]);
[test,timeboundaries{1}] = formattool(data(1,:), time, timelocks, binsize, maxLead, maxLag, method); % format the data and figure out the number of data points
DataRaster{1} = NaN([4,size(test,2),size(test,3)]);
DataRaster{1}(1,:,:) = test;
for iVar = 2:4
    DataRaster{1}(iVar,:,:) = formattool(data(iVar,:), time, timelocks, binsize, maxLead, maxLag, method);
end
[test,timeboundaries{2}] = formattool(data(5,:), time, timelocks, binsize, maxLead, maxLag, method); % process the second category
DataRaster{2} = NaN([2,size(test,2),size(test,3)]);
DataRaster{2}(1,:,:) = test;
DataRaster{2}(2,:,:) = formattool(data(6,:), time, timelocks, binsize, maxLead, maxLag, method);


% State the data using various numbers (2 to 7) of uniform width bins
MethodAssign = {1,1,'UniWB',{2};...
    1,2,'UniWB',{3};...
    1,3,'UniWB',{4};...
    1,4,'UniWB',{5};...
    2,1,'UniWB',{6};...
    2,2,'UniWB',{7}};
StatesRaster = data2states(DataRaster, MethodAssign);

% Perform the information calculation for all time bins
Method = 'Ent';
InfoResults = cell([2,1]);
for iDataCategory = 1:2
    nVar = size(StatesRaster{iDataCategory},1);
    nT = size(StatesRaster{iDataCategory},2);
    InfoResults{iDataCategory} = NaN([nVar,nT]);
    for iVar = 1:nVar
        for iT = 1:nT
            VariableIDs = {iDataCategory,iVar,iT};
            InfoResults{iDataCategory}(iVar,iT) = instinfo(StatesRaster, Method, VariableIDs);
        end
    end
end

figure
hold on
plot(mean(timeboundaries{1},2),InfoResults{1}(1,:),'r')
plot(mean(timeboundaries{1},2),InfoResults{1}(2,:),'g')
plot(mean(timeboundaries{1},2),InfoResults{1}(3,:),'k')
plot(mean(timeboundaries{1},2),InfoResults{1}(4,:),'b')
plot(mean(timeboundaries{2},2),InfoResults{2}(1,:),'c')
plot(mean(timeboundaries{2},2),InfoResults{2}(2,:),'m')
xlabel('Time')
ylabel('Information (bits)')
title('Entropy')
legend('Variable 1 (2 bins)','Variable 2 (3 bins)','Variable 3 (4 bins)','Variable 4 (5 bins)','Variable 5 (6 bins)','Variable 6 (7 bins)')