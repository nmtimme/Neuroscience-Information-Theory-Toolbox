%% Transfer Entropy Demo 1
% This demo goes through calculating the transfer entropy between 2
% continuous time series with time locked trials.

% Data Type: Continuous
% Number of Variables: 2
% Information Analysis: Transfer Entropy
% Trial Based or Single Trial: Trial Based
% Significance Testing: On
% Delays: Off

% We will generate data where the past state of variable 1 predicts the
% future state of variable 2 beyond what the past of variable 2 predicts
% about the future of variable 2
interaction = [zeros([1,10]),linspace(0,1,10),zeros([1,10])]; % the interaction will be off, than ramp up, then turn off again
interaction = repmat(interaction,[1,40]); % run the trial 40 times
copyLocs = rand(size(interaction)) < interaction; % find the times when the interaction made variable 2 copy the state of variable 1
data = randn([2,length(interaction)]); % generate the data
for it = 2:length(interaction) % impose self interaction from variable 2
    if rand < 0.1
        data(2,it) = data(2,it - 1);
    end
end
data(2,[false,copyLocs(1:(end - 1))]) = data(1,[copyLocs(1:(end - 1)),false]); % impose the copy states from variable 1 to variable 2 with a 1 time bin delay
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

% State the data using 4 uniform count bins
MethodAssign = {1,1,'UniCB',{4};...
    1,2,'UniCB',{4}};
StatesRaster = data2states(DataRaster, MethodAssign);

% Perform the information calculation for all time bins
Method = 'TE';
nT = size(StatesRaster,2);
InfoResults = NaN([nT,1]);
p = NaN([nT,1]);
for iT = 2:nT
    VariableIDs = {1,2,iT;... % Receiving variable in the future
        1,2,iT - 1;... % Receiving variable in the past
        1,1,iT - 1}; % Transmitting variable in the past
    [InfoResults(iT),p(iT)] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', 'on');
end

figure
plot(mean(timeboundaries,2),InfoResults)
xlabel('Time')
ylabel('Information (bits)')
title('Transfer Entropy')

figure
plot(mean(timeboundaries,2),p)
xlabel('Time')
ylabel('p-value')
set(gca,'YScale','log')
title('Transfer Entropy p-values')