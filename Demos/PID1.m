%% Partial Information Decomposition Demo 1
% This demo goes through calculating the partial information decomposition
% for time locked trials in 3 time series.

% Data Type: Discrete
% Number of Variables: 3
% Information Analysis: Partial Information Decomposition
% Trial Based or Single Trial: Trial Based
% Significance Testing: On
% Delays: Off

% For this demo, we will have an AND-gate operation that turns on and off
% in the trial
interaction = [zeros([1,10]),ones([1,10]),zeros([1,10])]; % the interaction will be off, than on, then turn off again
interaction = repmat(interaction,[1,40]); % run the trial 40 times
dataX1 = randi(2,[1,length(interaction)]); % random X1 variable
dataX2 = randi(2,[1,length(interaction)]); % random X2 variable
dataY = randi(2,[1,length(interaction)]);
dataY((interaction == 1) & (dataX1 == 2) & (dataX2 == 2)) = 2; % impose the AND operation
dataY((interaction == 1) & ~((dataX1 == 2) & (dataX2 == 2))) = 1; % impose the AND operation
time = 1:length(interaction); % assume the data are binned at 1 time unit
timelocks = 11:30:time(end); % lock the trials to the start of the interaction
binsize = 1; % leave the bins unchanged
maxLead = 10; % the interaction was off for 10 time steps prior to the time lock on each trial
maxLag = 20; % the interaction ramped up and then turned off for a total of 20 time steps after the time lock
method = 'ave'; % the method is not relevant because we did not change the bin size

% Format the data based on the trials
DataRaster = NaN([3,30,40]);
[DataRaster(1,:,:),timeboundaries] = formattool(dataY, time, timelocks, binsize, maxLead, maxLag, method);
DataRaster(2,:,:) = formattool(dataX1, time, timelocks, binsize, maxLead, maxLag, method);
DataRaster(3,:,:) = formattool(dataX2, time, timelocks, binsize, maxLead, maxLag, method);

% Assume the data are already in the appropriate states, so use the native
% stating
MethodAssign = []; % 
StatesRaster = data2states(DataRaster, MethodAssign);

% Perform the information calculation for all time bins
Method = '2PID';
nT = size(StatesRaster,2);
InfoResults = NaN([nT,4]);
p = NaN([nT,4]);
for iT = 1:nT
    VariableIDs = {1,1,iT;... % Y Variable
        1,2,iT;... % X1 Variable
        1,3,iT}; % X2 Variable
    [InfoResults(iT,:),p(iT,:)] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', 'on');
end

figure
hold on
plot(mean(timeboundaries,2),InfoResults(:,1),'Color','r')
plot(mean(timeboundaries,2),InfoResults(:,2),'Color','g')
plot(mean(timeboundaries,2),InfoResults(:,3),'Color','k')
plot(mean(timeboundaries,2),InfoResults(:,4),'Color','b')
legend('Redundancy','Unique X1','Unique X2','Synergy')
xlabel('Time')
ylabel('Information (bits)')
title('Partial Information Results')

figure
hold on
plot(mean(timeboundaries,2),p(:,1),'Color','r')
plot(mean(timeboundaries,2),p(:,2),'Color','g')
plot(mean(timeboundaries,2),p(:,3),'Color','k')
plot(mean(timeboundaries,2),p(:,4),'Color','b')
legend('Redundancy','Unique X1','Unique X2','Synergy')
xlabel('Time')
ylabel('p-value')
set(gca,'YScale','log')
title('Partial Information p-values')