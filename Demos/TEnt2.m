%% Transfer Entropy Demo 2
% This demo goes through calculating the transfer entropy between two
% discrete time series. In particular, this example works with neural spike
% trains.

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