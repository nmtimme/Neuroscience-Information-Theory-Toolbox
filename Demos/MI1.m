%% Mutual Information Demo 1
% This demo goes through calculating the mutual information for time locked
% trials from two discrete time series.

% Data Type: Discrete
% Number of Variables: 2
% Information Analysis: Mutual Information
% Trial Based or Single Trial: Trial Based
% Significance Testing: Off
% Delays: Off

% For this demo, we will generate discrete data (4 possible values) and
% vary the interaction between the variables. 
interaction = [zeros([1,10]),linspace(0,1,10),zeros([1,10])]; % the interaction will be off, than ramp up, then turn off again
interaction = repmat(interaction,[1,40]); % run the trial 40 times
copyLocs = rand(size(interaction)) < interaction; % find the times when the interaction made variable 2 copy the state of variable 1
data = randi(4,[2,length(interaction)]); % generate the data
data(2,copyLocs) = data(1,copyLocs); % impose the copy states
time = 1:length(interaction); % assume the data are binned at 1 time unit
timelocks = 11:30:time(end); % lock the trials to the start of the interaction strength ramp up
binsize = 1; % leave the bins unchanged
maxLead = 10; % the interaction was off for 10 time steps prior to the time lock on each trial
maxLag = 20; % the interaction ramped up and then turned off for a total of 20 time steps after the time lock
method = 'ave'; % the method is not relevant because we did not change the bin size

% Format the data based on the trials
DataRaster = NaN([2,30,40]);
[DataRaster(1,:,:),timeboundaries] = formattool(data(1,:), time, timelocks, binsize, maxLead, maxLag, method);
DataRaster(2,:,:) = formattool(data(2,:), time, timelocks, binsize, maxLead, maxLag, method);

% Assume the data are already in the appropriate states, so use the native
% stating
MethodAssign = []; % 
StatesRaster = data2states(DataRaster, MethodAssign);

% Perform the information calculation for all time bins
Method = 'PairMI';
nT = size(StatesRaster,2);
InfoResults = NaN([nT,1]);
for iT = 1:nT
    VariableIDs = {1,1,iT;1,2,iT};
    InfoResults(iT) = instinfo(StatesRaster, Method, VariableIDs);
end

figure
plot(mean(timeboundaries,2),InfoResults)
xlabel('Time')
ylabel('Information (bits)')
title('Mutual Information')