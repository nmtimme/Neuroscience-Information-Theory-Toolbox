%% demos - demonstrations of the analysis software
% Performs various information theoretic analyses using simple model data
% generated within. Each specific type of demo can be run individually.
% 
% Other m-files required: whole software package
% Subfunctions: none
% MAT-files required: 2PIDMats.mat, 3PIDMats.mat, TE3Redux.mat
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% May 2015; Last revision: 13-Jun-2017


%% Demo 1
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





%% Demo 2
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





%% Demo 3
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

        




%% Demo 4
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





%% Demo 5
% Data Type: Continuous
% Number of Variables: 2
% Information Analysis: Conditional Entropy
% Trial Based or Single Trial: Trial Based
% Significance Testing: N/A
% Delays: Off


% For this demo, we will generate continuous data and vary the interaction 
% between the variables. 
interaction = [zeros([1,10]),linspace(0,1,10),zeros([1,10])]; % the interaction will be off, than ramp up, then turn off again
interaction = repmat(interaction,[1,40]); % run the trial 40 times
copyLocs = rand(size(interaction)) < interaction; % find the times when the interaction made variable 2 copy the state of variable 1
data = randn([2,length(interaction)]); % generate the data
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

% State the data using 4 uniform width bins
MethodAssign = {1,1,'UniWB',{4};...
    1,2,'UniWB',{4}};
StatesRaster = data2states(DataRaster, MethodAssign);

% Perform the conditional entropy of variable 1 conditioned on variable 2
% at all time bins
Method = 'CondEnt1';
nT = size(StatesRaster,2);
InfoResults = NaN([nT,1]);
for iT = 1:nT
    VariableIDs = {1,1,iT;...
        1,2,iT};
    InfoResults(iT) = instinfo(StatesRaster, Method, VariableIDs);
end

figure
plot(mean(timeboundaries,2),InfoResults)
xlabel('Time')
ylabel('Information (bits)')
title('Entropy')


%% Demo 6
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





%% Demo 7
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





%% Demo 8
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





%% Demo 9
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





%% Demo 10
% Data Type: Discrete (Neuron Action Potentials)
% Number of Variables: 2
% Information Analysis: Transfer Entropy
% Trial Based or Single Trial: Single Trial
% Significance Testing: On
% Delays: On

% For this demo, make a driving neuron with random Poisson spikes. The
% receiving neuron will have an identical background firing rate, but 15
% time bins after the transmitter spikes, the receiver will have an
% increased likelihood to spike.
nData = 10000;
dataTrans = zeros([1,nData]);
dataTrans(rand([1,nData]) < 0.1) = 1;
pRec = 0.1*ones([1,nData]);
pRec(16:end) = pRec(16:end) + 0.4*dataTrans(1:(end - 15));
dataRec = zeros([1,nData]);
dataRec(rand([1,nData]) < pRec) = 1;

% Format the data to increase the bin size
time = 1:nData; % assume the data are binned at 1 time unit
timelocks = 1; % There is only one trial
binsize = 5; % rebin the data
maxLead = 0; % since the timelock is the first data point, there are no data to consider earlier
maxLag = nData; % we want to include the whole recording
method = 'count'; % we will count the number of spikes in the bins
DataRaster = NaN([2,nData/binsize]);
[DataRaster(1,:),timeboundaries] = formattool(dataTrans, time, timelocks, binsize, maxLead, maxLag, method);
DataRaster(2,:) = formattool(dataRec, time, timelocks, binsize, maxLead, maxLag, method);

% State the data using 2 uniform count bins
MethodAssign = {1,1,'UniCB',{2};...
    1,2,'UniCB',{2}};
StatesRaster = data2states(DataRaster, MethodAssign);

% Perform the information calculation with the wrong delay
Method = 'TE';
VariableIDs = {1,2,2;... % Receiving variable in the future
    1,2,1;... % Receiving variable in the past
    1,1,1}; % Transmitting variable in the past
[InfoResults,p] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', 'on');

disp(['Incorrect Delay Transfer Entropy = ',num2str(InfoResults),' bits (p = ',num2str(p(1)),')'])

% Perform the information calculation with the correct delay
Method = 'TE';
VariableIDs = {1,2,4;... % Receiving variable in the future
    1,2,3;... % Receiving variable in the past
    1,1,1}; % Transmitting variable in the past
[InfoResults,p] = instinfo(StatesRaster, Method, VariableIDs, 'MCOpt', 'on');

disp(['Correct Delay Transfer Entropy = ',num2str(InfoResults),' bits (p = ',num2str(p(1)),')'])






%% Demo 11
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


%% Demo 12
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

