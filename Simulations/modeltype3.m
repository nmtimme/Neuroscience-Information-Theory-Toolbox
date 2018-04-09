%% modeltype3 - full Izhikevich network
% Generates data for the full Izhikevich neuron network. This network is
% precisely the 1000 neuron network described in Izhikevich 2003. In
% addition to the spontaneous activity produced by the network, square
% current pulses are supplied to subsets of excitatory neurons. The neurons 
% in this model do not have membrane noise.
%
% E. M. Izhikevich, Simple model of spiking neurons, IEEE Transactions on
% Neural Networks 14 (6), 2003. 
%
%
% Syntax: [spkTimes,neuronID,stim,conmat,modelParams,neuronPos] = modeltype3(stimType,varargin)
%
% Inputs:
%   stimType - determines the type of stimulation to be applied (string).
%       Options are:
%           'OneStim': A single square stimulus pulse is applied to
%             nStimRec randomly selected excitatory neurons.
%           'TwoStim': Square stimulus pulses are applied to two sets of
%             randomly selected excitatory neurons. Every combination of on
%             and off for both stimuli is applied. 
% 
% Optional Inputs:
%   modeltype3(..., 'nStimRec', nStimRec) - sets the number of excitatory
%       neurons in each group that receives stimulus current (double,
%       default = 40). The neurons in each group are randomly selected from
%       all excitatory neurons, so there can be overlap between the groups.
%   modeltype3(..., 'pulseSize', pulseSize) - sets the size of the square
%       current pulses. The current is unitless in this model (double,
%       default = 50).
%   modeltype3(..., 'interstimT', interstimT) - sets the time between
%       stimulation pulses in milliseconds (double, default = 1000).
%   modeltype3(..., 'intrastimT', intrastimT) - sets the duration of the
%       stimulation pulses in milliseconds (double, default = 100). 
%   modeltype3(..., 'nStimSets', nStimSets) - the number of stimulation
%       sets to apply (integer, default = 50).
%   modeltype3(..., 'pulseValueSTD', pulseValueSTD) - sets the standard
%       deviation of the Gaussian white noise added to the input current
%       received by the model neurons (double, default = 0.05). 
%   modeltype3(..., 'pulseOnsetSTD', pulseOnsetSTD) - sets the standard
%       deviation of the Gaussian white noise added to the onset and offset
%       time of the stimuli pulses received by the model neurons in
%       milliseconds (double, default = 5). 
%   modeltype3(..., 'stabilizeTime', stabilizeTime) - sets the time in
%       milliseconds for the network to run without stimuli to stabilize
%       (double, default = 3000). 
%   modeltype3(..., 'conMode', conMode') - sets the connectivity (string,
%       default = 'AllToAll'). Options:
%           'AllToAll': All neurons are connected to all other neurons with
%               random uniform weights per Izhikevich's original network. 
%               If this option is selected, neuronPos = [].
%           'Spatial': The neurons are placed on a 2-D surface with
%               periodic boundary conditions and their connectivity
%               probability is based on the shortest distance between the
%               neurons. The likelihood to be connected decays
%               exponentially with distance. The synaptic weight is
%               increased to produce the same average total synaptic weight
%               as would be found in the original all-to-all connectivity.
%               The 2-D surface is a square with side length 1.
%   modeltype3(..., 'spatialDecay', spatialDecay) - sets the exponential
%       decay in the spatial network connectivity (double, default = 0.05).
%       This parameter represents the Cartesian distance for which the
%       likelihood for two neurons to be connected is 0.5. If conMode =
%       'AllToAll', then spatialDecay has no effect. 
%   modeltype3(..., 'weightScaling', weightScaling) - sets the weight
%       rescaling for the spatial network connectivity (double, default =
%       0.5). weightScaling represents the ratio of the total synaptic
%       weight in the spatial network to the total synaptic weight in the
%       original all-to-all network. If conMode = 'AllToAll', then
%       weightScaling has no effect.
%   modeltype3(..., 'spatialDist', spatialDist) - sets the distribution of
%       the stimulation points in the spatial connection mode (string, 
%       default - 'Uniform'). If conMode = 'AllToAll', then spatialDist has
%       no effect. Options:
%           'Uniform': the stimulation points are randomly distributed
%               throughout the 1 by 1 2-D surface. 
%           'Point': the stimulation points are concentrated near specific
%               point(s) in the 1 by 1 2-D surface. When stimType =
%               'OneStim', the point is (0.5,0.5). When stimType =
%               'TwoStim', the points are (0.25,0.25) for stim 1 and
%               (0.75,0.75) for stim 2.
%           'Line': the stimulation points are concentrated near a line in
%               1 by 1 2-D surface. When stimType = 'OneStim', the line is
%               x = 0.5. When stimType = 'TwoStim', the lines are x = 0.25
%               for stim 1 and x = 0.75 for stim 2.
%
% Outputs:
%   spkTimes (1000 by 1 cell array, each element is a double vector) - list 
%       of times when a neuron spiked in milliseconds. 
%   neuronID (1000 by 1 cell array, each element is a string) - list of
%       neuron types for the correspond element of spkTimes. Options
%       include:
%           '1': The neuron received inputs from stimulus 1.
%           '2': The neuron received inputs from stimulus 2. (Only applies
%               to stimType = 'TwoStim')
%           '12': The neuron received inputs from stimuli 1 and 2. (Only
%               applies to stimType = 'TwoStim')
%           'Ex': The neuron was an excitatory neuron that did not receive
%               direct stimulus input.
%           'In': The neuron was an inhibitory neuron.
%   stim (structure) - a list of various quantities associated with the
%       stimuli:
%       stim.times (double vector) - the onset times for the stimuli in
%           milliseconds. Note: these times are recorded prior to adding
%           noise to the input current or the application of other
%           manipulations.
%       stim.types (integer vector) - the type of stimuli that was applied
%           at the time listed in stim.times. 
%       stim.currentOriginal (2 by number of time steps double vector) - 
%           the original clean current for stimulus 1 (top row) and 
%           stimulus 2 (bottom row).
%       stim.currentRec(2 by number of time steps double vector) - the
%           actual stimulus current after adding noise for stimulus 1 (top
%           row) and stimulus 2 (bottomw row).
%   conmat (1000 by 1000 double array) - connectivity matrix. conmat(i,j)
%       represents the weight from neuron i to neuron j.
%   modelParams (structure) - a list of the parameters used in the model.
%       Fields are named using the variable input names listed above.
%
% Examples:
%   Default Stimulated Izhikevich 1000 Neuron Model
%       [spkTimes,stim,modelParams] = modeltype3;
%   Default Stimulated Izhikevich 1000 Neuron Model with 100 stimulation
%   sets instead of 50.
%       [spkTimes,stim,modelParams] = modeltype3('nStimSets',100);
%   Default Stimulated Izhikevich 1000 Neuron Model with 10 seconds of
%   stabilization time instead of 2 seconds.
%       [spkTimes,stim,modelParams] = modeltype3('stabilizeTime',10000);
%   
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% March 2016; Last revision: 2-Mar-2016


function [spkTimes,neuronID,stim,conmat,modelParams,neuronPos] = modeltype3(stimType,varargin)
%% Parse command line for parameters

% Set default values
nStimRec = 40;
pulseSize = 50;
interstimT = 1000;
intrastimT = 10;
nStimSets = 50;
pulseValueSTD = 0.05;
pulseOnsetSTD = 1;
stabilizeTime = 3000;
neuronPos = [];
conMode = 'AllToAll';
spatialDecay = 0.05;
weightScaling = 0.5;
spatialDist = 'Uniform';

% Loop through the inputs
iVarArg = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'nStimRec',            nStimRec = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'pulseSize',           pulseSize = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'interstimT',          interstimT = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'intrastimT',          intrastimT = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'nStimSets',           nStimSets = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'pulseValueSTD',       pulseValueSTD = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'pulseOnsetSTD',       pulseOnsetSTD = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'stabilizeTime',       stabilizeTime = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'conMode',             conMode = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'spatialDecay',        spatialDecay = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'weightScaling',       weightScaling = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'spatialDist',         spatialDist = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(MODELTYPE3) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Set the parameters from the original Izhikevich model

% Excitatory neurons    Inhibitory neurons
Ne=800;                 Ni=200;
re=rand(Ne,1);          ri=rand(Ni,1);
a=[0.02*ones(Ne,1);     0.02+0.08*ri];
b=[0.2*ones(Ne,1);      0.25-0.05*ri];
c=[-65+15*re.^2;        -65*ones(Ni,1)];
d=[8-6*re.^2;           2*ones(Ni,1)];
S=[0.5*rand(Ne+Ni,Ne),  -rand(Ne+Ni,Ni)];

% Error check the number of neurons
if (Ne + Ni) ~= 1000
    error('We were expecting 1000 total neurons.')
end


%% Make the Input Current, Add Noise, and Make the Membrane Noise

if strcmp(stimType,'TwoStim')
    
    % Reorganize the pulse size information
    pulseSize = [0,pulseSize,0,pulseSize;0,0,pulseSize,pulseSize];
    
    % Find the number of pulses in one stimulation set
    nPulses = size(pulseSize,2);
    
    % Find the total length of the recording in milliseconds
    T = nStimSets*nPulses*intrastimT + nStimSets*nPulses*interstimT + interstimT + stabilizeTime;
    
    % Make the stimulus currents
    I1 = [];
    I2 = [];
    stims1 = [];
    stims2 = [];
    for iPulse = 1:nPulses
        I1 = [I1,zeros([1,interstimT]),pulseSize(1,iPulse)*ones([1,intrastimT])];
        stims1 = [stims1,NaN([1,interstimT]),pulseSize(1,iPulse),NaN([1,intrastimT - 1])];
        I2 = [I2,zeros([1,interstimT]),pulseSize(2,iPulse)*ones([1,intrastimT])];
        stims2 = [stims2,NaN([1,interstimT]),pulseSize(2,iPulse),NaN([1,intrastimT - 1])];
    end
    I1 = repmat(I1,[1,nStimSets]);
    I1 = [zeros([1,stabilizeTime]),I1,zeros([1,interstimT])];
    I2 = repmat(I2,[1,nStimSets]);
    I2 = [zeros([1,stabilizeTime]),I2,zeros([1,interstimT])];
    
    % Error check stimulation length
    if length(I1) ~= T
        error('Stimulation length error.')
    elseif length(I2) ~= T
        error('Stimulation length error.')
    end
    
    stims1 = repmat(stims1,[1,nStimSets]);
    stims1 = [NaN([1,stabilizeTime]),stims1];
    stims2 = repmat(stims2,[1,nStimSets]);
    stims2 = [NaN([1,stabilizeTime]),stims2];
    
    [B,I,J] = unique(stims1(isfinite(stims1)));
    stims1(isfinite(stims1)) = J;
    [B,I,J] = unique(stims2(isfinite(stims2)));
    stims2(isfinite(stims2)) = J;
    
    stim = struct;
    stim.times = cell([2,1]);
    stim.times{1} = find(isfinite(stims1));
    stim.times{2} = find(isfinite(stims2));
    stim.types = cell([2,1]);
    stim.types{1} = stims1(isfinite(stims1));
    stim.types{2} = stims2(isfinite(stims2));
    
    % Keep a copy of the original stimulation current
    stim.currentOriginal = [I1;I2];
    
    % Add noise to the pulse onset and offset
    transpoints = find(diff(I1) ~= 0);
    noise = round(pulseOnsetSTD*randn(size(transpoints)));
    for itrans = 1:length(transpoints)
        if noise(itrans) > 0
            I1((transpoints(itrans) + 1):(transpoints(itrans) + noise(itrans))) = I1(transpoints(itrans));
        elseif noise(itrans) < 0
            I1((transpoints(itrans) + noise(itrans) + 1):transpoints(itrans)) = I1(transpoints(itrans) + 1);
        end
    end
    transpoints = find(diff(I2) ~= 0);
    noise = round(pulseOnsetSTD*randn(size(transpoints)));
    for itrans = 1:length(transpoints)
        if noise(itrans) > 0
            I2((transpoints(itrans) + 1):(transpoints(itrans) + noise(itrans))) = I2(transpoints(itrans));
        elseif noise(itrans) < 0
            I2((transpoints(itrans) + noise(itrans) + 1):transpoints(itrans)) = I2(transpoints(itrans) + 1);
        end
    end
    
    % Add noise to the input current
    I1 = I1 + pulseValueSTD*randn(size(I1));
    I2 = I2 + pulseValueSTD*randn(size(I2));
    
    % Keep a copy of the stimulation current received by the neuron
    stim.currentRec = [I1;I2];

elseif strcmp(stimType,'OneStim')
    
    % Reorganize the pulse size information
    pulseSize = [0,pulseSize];
    
    % Find the number of pulses in one stimulation set
    nPulses = size(pulseSize,2);
    
    % Find the total length of the recording in milliseconds
    T = nStimSets*nPulses*intrastimT + nStimSets*nPulses*interstimT + interstimT + stabilizeTime;
    
    % Make the stimulus currents
    I1 = [];
    stims1 = [];
    for iPulse = 1:nPulses
        I1 = [I1,zeros([1,interstimT]),pulseSize(1,iPulse)*ones([1,intrastimT])];
        stims1 = [stims1,NaN([1,interstimT]),pulseSize(1,iPulse),NaN([1,intrastimT - 1])];
    end
    I1 = repmat(I1,[1,nStimSets]);
    I1 = [zeros([1,stabilizeTime]),I1,zeros([1,interstimT])];
    
    % Error check stimulation length
    if length(I1) ~= T
        error('Stimulation length error.')
    end
    
    stims1 = repmat(stims1,[1,nStimSets]);
    stims1 = [NaN([1,stabilizeTime]),stims1];
    
    [B,I,J] = unique(stims1(isfinite(stims1)));
    stims1(isfinite(stims1)) = J;
    
    stim = struct;
    stim.times = find(isfinite(stims1));
    stim.types = stims1(isfinite(stims1));

    
    % Keep a copy of the original stimulation current
    stim.currentOriginal = [I1];
    
    % Add noise to the pulse onset and offset
    transpoints = find(diff(I1) ~= 0);
    noise = round(pulseOnsetSTD*randn(size(transpoints)));
    for itrans = 1:length(transpoints)
        if noise(itrans) > 0
            I1((transpoints(itrans) + 1):(transpoints(itrans) + noise(itrans))) = I1(transpoints(itrans));
        elseif noise(itrans) < 0
            I1((transpoints(itrans) + noise(itrans) + 1):transpoints(itrans)) = I1(transpoints(itrans) + 1);
        end
    end
    
    % Add noise to the input current
    I1 = I1 + pulseValueSTD*randn(size(I1));
    
    % Keep a copy of the stimulation current received by the neuron
    stim.currentRec = I1;
    
    % Make the necessary empty vector for stimulus two
    I2 = zeros(size(I1));
    In2 = [];

end


%% If Necessary, Make the Connectivity Spatial Dependent

if strcmp(conMode,'Spatial')
    
    % Randomly select neuron positions
    neuronPos = rand([1000,2]);
    
    % Make the mirror positions for periodic boundary conditions
    mirrorPos = zeros([1000,2,9]);
    iSurface = 1;
    for i = [-1,0,1]
        for j = [-1,0,1]
            mirrorPos(:,1,iSurface) = neuronPos(:,1) + i;
            mirrorPos(:,2,iSurface) = neuronPos(:,2) + j;
            iSurface = iSurface + 1;
        end
    end
    
    % Calculate the minimum distance between all neuron pairs
    distanceMat = zeros([1000,1000]);
    for i = 1:999
        for j = (i + 1):1000
            distanceMat(i,j) = min(sqrt(sum((mirrorPos(j,:,:) - neuronPos(i,:,ones([1,9]))).^2,2)));
            distanceMat(j,i) = distanceMat(i,j);
        end
    end
    
    % Calculate the decay exponent
    alpha = -log(0.5)/spatialDecay;
    
    % Convert the distance matrix to a probability matrix
    probMat = exp(-alpha.*distanceMat);
    
    % Create the new connectivity matrix
    Snew = S;
    Snew(rand([1000,1000]) > probMat) = 0;
    
    % Rescale the weights to match the same total from the original
    % all-to-all matrix
    totalExWeightOld = sum(sum(S(:,1:Ne),1),2);
    totalExWeightNew = sum(sum(Snew(:,1:Ne),1),2);
    totalInWeightOld = sum(sum(S(:,(Ne + 1):1000),1),2);
    totalInWeightNew = sum(sum(Snew(:,(Ne + 1):1000),1),2);
    Snew(:,1:Ne) = weightScaling*(totalExWeightOld/totalExWeightNew).*Snew(:,1:Ne);
    Snew(:,(Ne + 1):1000) = weightScaling*(totalInWeightOld/totalInWeightNew).*Snew(:,(Ne + 1):1000);
    
    % Relabel the connection matrix
    S = Snew;
    
    
end


%% Assign the Stimulation Neurons


if strcmp(stimType,'TwoStim')
    
    if strcmp(conMode,'Spatial') && strcmp(spatialDist,'Point')
        
        % Find the excitatory neurons closest to (0.25,0.25)
        distance = sqrt(sum((neuronPos(1:Ne,:) - repmat([0.25,0.25],[Ne,1])).^2,2));
        [waste,I] = sort(distance,1,'ascend');
        In1 = I(1:nStimRec);
        
        % Find the excitatory neurons closest to (0.75,0.75)
        distance = sqrt(sum((neuronPos(1:Ne,:) - repmat([0.75,0.75],[Ne,1])).^2,2));
        [waste,I] = sort(distance,1,'ascend');
        In2 = I(1:nStimRec);
        
    elseif strcmp(conMode,'Spatial') && strcmp(spatialDist,'Line')
        
        % Find the excitatory neurons closest to x = 0.25
        distance = abs(neuronPos(1:Ne,1) - 0.25);
        [waste,I] = sort(distance,1,'ascend');
        In1 = I(1:nStimRec);
        
        % Find the excitatory neurons closest to x = 0.75
        distance = abs(neuronPos(1:Ne,1) - 0.75);
        [waste,I] = sort(distance,1,'ascend');
        In2 = I(1:nStimRec);
        
    else
        
        % Randomly select stimulus input sets
        tempRand = randperm(Ne);
        In1 = tempRand(1:nStimRec);
        tempRand = randperm(Ne);
        In2 = tempRand(1:nStimRec);
        
    end
    
    % Make the list of neuron identifiers
    neuronID = cell([1000,1]);
    for iNeuron = 1:1000
        if ismember(iNeuron,In1)
            if ismember(iNeuron,In2)
                neuronID{iNeuron} = '12';
            else
                neuronID{iNeuron} = '1';
            end
        elseif ismember(iNeuron,In2)
            neuronID{iNeuron} = '2';
        elseif iNeuron <= 800
            neuronID{iNeuron} = 'Ex';
        else
            neuronID{iNeuron} = 'In';
        end
    end

elseif strcmp(stimType,'OneStim')
    
    if strcmp(conMode,'Spatial') && strcmp(spatialDist,'Point')
        
        % Find the excitatory neurons closest to (0.5,0.5)
        distance = sqrt(sum((neuronPos(1:Ne,:) - repmat([0.5,0.5],[Ne,1])).^2,2));
        [waste,I] = sort(distance,1,'ascend');
        In1 = I(1:nStimRec);
        
    elseif strcmp(conMode,'Spatial') && strcmp(spatialDist,'Line')
        
        % Find the excitatory neurons closest to x = 0.5
        distance = abs(neuronPos(1:Ne,1) - 0.5);
        [waste,I] = sort(distance,1,'ascend');
        In1 = I(1:nStimRec);
        
    else
        
        % Randomly select stimulus input sets
        tempRand = randperm(Ne);
        In1 = tempRand(1:nStimRec);
        
    end
    
    % Make the list of neuron identifiers
    neuronID = cell([1000,1]);
    for iNeuron = 1:1000
        if ismember(iNeuron,In1)
            neuronID{iNeuron} = '1';
        elseif iNeuron <= 800
            neuronID{iNeuron} = 'Ex';
        else
            neuronID{iNeuron} = 'In';
        end
    end

end



%% Run the Model

% Preallocate space to record the spikes
spkTimes = cell([1000,1]);
for i = 1:1000
    spkTimes{i} = NaN([1,1000]);
end
spkInd = ones([1000,1]);

v=-65*ones(Ne+Ni,1);    % Initial values of v
u=b.*v;                 % Initial values of u
spkTimes = cell([1000,1]);

for t=1:T            % simulation of 1000 ms
    I=[5*randn(Ne,1);2*randn(Ni,1)]; % thalamic input
    fired=find(v>=30);    % indices of spikes
    
    % Record the spikes
    for iSpike = 1:length(fired)
        if spkInd(fired(iSpike)) > length(spkTimes{fired(iSpike)})
            spkTimes{fired(iSpike)} = [spkTimes{fired(iSpike)},NaN([1,1000])];
        end
        spkTimes{fired(iSpike)}(spkInd(fired(iSpike))) = t;
        spkInd(fired(iSpike)) = spkInd(fired(iSpike)) + 1;
    end
    
    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);
    I=I+sum(S(:,fired),2);
    
    % Add the stimuli
    I(In1) = I(In1) + I1(t);
    I(In2) = I(In2) + I2(t);
    
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u);                 % stability
end;

% Clean up the spike times
for iNeuron = 1:1000
    spkTimes{iNeuron}(~isfinite(spkTimes{iNeuron})) = [];
end



%% Organize the Results

% Take the transpose of the S matrix to get the connectivity matrix
conmat = S';

% Organize the parameters for this model
modelParams = struct;
modelParams.interstimT = interstimT;
modelParams.intrastimT = intrastimT;
modelParams.nStimSets = nStimSets;
modelParams.pulseValueSTD = pulseValueSTD;
modelParams.pulseOnsetSTD = pulseOnsetSTD;
modelParams.pulseSize = pulseSize;
modelParams.stabilizeTime = stabilizeTime;
modelParams.nStimRec = nStimRec;
modelParams.conMode = conMode;
modelParams.spatialDecay = spatialDecay;
modelParams.weightScaling = weightScaling;
modelParams.spatialDist = spatialDist;


end