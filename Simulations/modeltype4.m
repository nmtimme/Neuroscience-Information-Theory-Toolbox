%% modeltype4 - 4 neuron logic gate-like Izhikevich neuron system
% Generates data for a small Izhikevich neuron network designed to
% reproduce basic logic gate-like operations. A pulse stimulis is applied 
% to two neurons. Those neurons send connections to a third and fourth 
% neuron. The third neuron also sends a connection to the fourth neuron.
% The fourth neuron also receives a background excitatory or inhibitory
% constant current (approximating activity in the rest of the network). The
% types of neurons, the connection weights, and the background constant
% current can be varied to produce various logic gate-like operations. The
% fourth neurons serves as the output.
%
% Syntax: [mV,spkTimes,memNoise,current,stim,modelParams] = modeltype4(pulseSize, weights, neuronTypes, backgroundCurrent, varargin)
%
% Inputs:
%   pulseSize - sets the sizes in pA of the stimulus pulse square waves
%       (double). The input neurons will receive every combination of
%       stimulus pulses ((0,0), (pulseSize,0), (0, pulseSize),
%       (pulseSize,pulseSize)).
%   weights - sets the maximum current in pA supplied from presynaptic
%       spikes to postsynaptic neurons (1 by 5 double array). Weight
%       assignments are as follows:
%           weights(1): first stimulus to output neuron
%           weights(2): second stimulus to output neuron
%           weights(3): first stimulus to third neuron
%           weights(4): second stimuls to third neuron
%           weights(5): third neuron to output neuron
%       Note: if inhibitory neurons are selected for either stimulus neuron, 
%       the corresponding weight must be negative (see neuronTypes).
%   neuronTypes - sets the types of the neurons (1 by 4 cellular array with
%       each element being a string). Options: 'Ex' for excitatory RS
%       neurons, 'Inh' for inhibitory FS neurons. neuronTypes{1} applies to
%       first stimulus neuron, neuronTypes{2} applies to the second
%       stimulus neuron, neuronTypes{3} applies to the third neuron,
%       neuronTypes{4} applies to the output neuron.
%   backgroundCurrent - sets the strength in pA of the background current
%       applied to the output neuron (double). This can be positive for
%       background excitatory input or negative for background inhibitory
%       input.
% 
% Optional Inputs:
%   modeltype4(..., 'noiseScale', noiseScale) - sets the membrane potential
%       noise strength. noiseScale (1 by 4 double vector) sets the strength 
%       (default = 0.1 for all neurons) and is the ratio of standard 
%       deviations of the noise signal.  The noise is modeled as 1/f noise 
%       bound below at 10 Hz with constant power noise from 1 Hz to 10 Hz. 
%   modeltype4(..., 'interstimT', interstimT) - sets the time between
%       stimulation pulses in milliseconds (double, default = 500).
%   modeltype4(..., 'intrastimT', intrastimT) - sets the duration of the
%       stimulation pulses in milliseconds (double, default = 500). 
%   modeltype4(..., 'nStimSets', nStimSets) - the number of stimulation
%       sets to apply (integer, default = 50).
%   modeltype4(..., 'tau', tau) - sets the bin size or sampling in
%       milliseconds (double, default = 0.1).
%   modeltype4(..., 'pulseValueSTD', pulseValueSTD) - sets the standard
%       deviation of the Gaussian white noise added to the input current
%       received by the model neuron in pA (double, default = 5). 
%   modeltype4(..., 'pulseOnsetSTD', pulseOnsetSTD) - sets the standard
%       deviation of the Gaussian white noise added to the onset and offset
%       time of the stimuli pulses received by the model neuron in
%       milliseconds (double, default = 5). 
%   modeltype4(..., 'gammaMean', gammaMean) - sets the mean of the gamma
%       pdf used to model current injection to postsynaptic neurons from
%       presynaptic spikes (1 by 4 double vector, default = 30*ones([1,4])). 
%       gammaMean(iNeuron) applies to all connections originating from
%       neuron iNeuron.
%   modeltype4(..., 'gammaSTD', gammaSTD) - sets the standard deviation of
%       the gamma pdf used to model current injection to postsynaptic
%       neurons from presynaptic spikes (1 by 4 double vector, default =
%       20*ones([1,4])). gammaSTD(iNeuron) applies to all connections
%       originating from neuron iNeuron.
%   modeltype4(..., 'varStim', varStim) - sets the probability for each
%       stimulation pattern (1 by 4 double vector, default =
%       0.25*ones([1,4])). Correspondence between varStim and the
%       stimulation pattern is given below:
%           (stim A = 0, stim B = 0): varStim(1)
%           (stim A = pulseSize, stim B = 0): varStim(2)
%           (stim A = 0, stim B = pulseSize): varStim(3)
%           (stim A = pulseSize, stim B = pulseSize): varStim(4)
%       Note: the precise number of presentations of a given stimulation
%       pattern will be round(varStim(i)*nStimSets*4).
%
% Outputs:
%   Note: for all output variables containing information about individual
%       neurons, the first element corresponds to stimulus neuron 1, the
%       second element corresponds to stimulus neuron 2, and the third
%       element corresponds to the third neuron, and the fourth element
%       corresponds to the output neuron.
%   mV (4 by number of time steps double array) - membrane potential of the 
%       model neurons in mV. 
%   spkTimes (4 by 1 cell array, each element is a double vector) - list of 
%       times when a neuron spiked in milliseconds. 
%   memNoise (4 by number of time steps double vector) - the membrane noise 
%       of the model neurons in mV. 
%   current (4 by number of time steps double vector) - the input current
%       to each neuron in pA. 
%   stim (structure) - a list of various quantities associated with the
%       stimuli:
%       stim.times (double vector) - the onset times for the stimuli in
%           milliseconds. Note: these times are recorded prior to adding
%           noise to the input current or the application of other
%           manipulations.
%       stim.types (integer vector) - the type of stimuli that was applied
%           at the time listed in stim.times. 
%       stim.currentOriginal (double vector) - the original clean current.
%       stim.currentRec (double vector) - the current after adding noise
%           and other possible manipulations.
%   modelParams (structure) - a list of the parameters used in the model.
%       Fields are named using the variable input names listed above.
%
% Examples:
%   
% Other m-files required: inverseFNoise, simpleModel_RS, simpleModel_FSI
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% March 2016; Last revision: 22-Mar-2016


function [mV,spkTimes,memNoise,current,stim,modelParams] = modeltype4(pulseSize, weights, neuronTypes, backgroundCurrent, varargin)
%% Parse command line for parameters


% Set default values
noiseScale = 0.1*ones([1,4]);
interstimT = 500;
intrastimT = 500;
nStimSets = 50;
tau = 0.1;
pulseValueSTD = 5;
pulseOnsetSTD = 5;
gammaMean = 30*ones([1,4]);
gammaSTD = 20*ones([1,4]);
varStim = 0.25*ones([1,4]);




iVarArg = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'noiseScale',          noiseScale = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'interstimT',          interstimT = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'intrastimT',          intrastimT = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'nStimSets',           nStimSets = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'tau',                 tau = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'pulseValueSTD',       pulseValueSTD = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'pulseOnsetSTD',       pulseOnsetSTD = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'gammaMean',           gammaMean = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'gammaSTD',            gammaSTD = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'varStim',             varStim = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(MODELTYPE4) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Do Some Preprocessing

% Error check the inhibitory weight
if strcmp(neuronTypes{1},'Inh') && ((weights(1) > 0) || (weights(3) > 0))
    error('The inhibitory weights from stimulus 1 neuron should be less than or equal to zero.')
end
if strcmp(neuronTypes{2},'Inh') && ((weights(2) > 0) || (weights(4) > 0))
    error('The inhibitory weights from stimulus 2 neuron should be less than or equal to zero.')
end
if strcmp(neuronTypes{3},'Inh') && (weights(5) > 0)
    error('The inhibitory weights from neuron three should be less than or equal to zero.')
end

% Error check that the stimulus pattern variation sums to 1
if sum(varStim) ~= 1
    error('The variable stimulus pattern must sum to 1.')
end

% Reorganize the pulse size information
pulseSize = [0,pulseSize,0,pulseSize;0,0,pulseSize,pulseSize];

% Find the number of pulses in one stimulation set
nPulses = size(pulseSize,2);

% Find the total length of the recording in milliseconds
T = nStimSets*nPulses*intrastimT + nStimSets*nPulses*interstimT + interstimT;

% Find the total number of time bins
n = round(T/tau);

%% Make the Input Current, Add Noise, and Make the Membrane Noise

% Identify and order the stimuli
nStims = 4*nStimSets;
jointStims = zeros([4,nStims]);
for iStim = 1:4
    jointStims(iStim,1:round(nStims*varStim(iStim))) = iStim;
end
jointStims = jointStims(:);
jointStims(jointStims == 0) = [];

% Figure out the stimulus times
stimTimes = (interstimT + tau):(interstimT + intrastimT):T;
stim = struct;
stim.times = cell([2,1]);
stim.times{1} = stimTimes;
stim.times{2} = stimTimes;

% Make the stimulus currents
I = cell([4,1]);
for iNeuron = 1:4
    I{iNeuron} = zeros([1,n]);
end
time = (1:n)*tau;
for iCur = 1:2
    for iStim = 1:nStims
        I{iCur}((time >= stim.times{iCur}(iStim)) & (time < (stim.times{iCur}(iStim) + intrastimT))) = pulseSize(iCur,jointStims(iStim));
    end
end

% Record the individual stimulus types
stim.types = cell([2,1]);
stim.types{1} = jointStims;
stim.types{1}(stim.types{1} == 3) = 1;
stim.types{1}(stim.types{1} == 4) = 2;
stim.types{2} = jointStims;
stim.types{2}(stim.types{2} == 2) = 1;
stim.types{2}(stim.types{2} == 3) = 2;
stim.types{2}(stim.types{2} == 4) = 2;

% Keep a copy of the original stimulation current
stim.currentOriginal = [I{1};I{2}];

% Add noise to the pulse onset and offset
for iCur = 1:2
    transpoints = find(diff(I{iCur}) ~= 0);
    noise = round((pulseOnsetSTD/tau)*randn(size(transpoints)));
    for itrans = 1:length(transpoints)
        if noise(itrans) > 0
            I{iCur}((transpoints(itrans) + 1):(transpoints(itrans) + noise(itrans))) = I{iCur}(transpoints(itrans));
        elseif noise(itrans) < 0
            I{iCur}((transpoints(itrans) + noise(itrans) + 1):transpoints(itrans)) = I{iCur}(transpoints(itrans) + 1);
        end
    end
end

% Add noise to the input current
for iCur = 1:2
    I{iCur} = I{iCur} + pulseValueSTD*randn(size(I{iCur}));
end

% Keep a copy of the stimulation current received by the neuron
stim.currentRec = [I{1};I{2}];

% Make the membrane noise
memNoise = NaN([4,n]);
for i = 1:4
    memNoise(i,:) = inverseFNoise(n, 1, 1000/tau, 'minF', 10, 'constLow', 1);
    memNoise(i,:) = memNoise(i,:)./std(memNoise(i,:));
    memNoise(i,:) = memNoise(i,:)*noiseScale(i);
end

% Make the necessary gamma pdfs for the current injected by spikes
gammaMean = gammaMean./tau;
gammaSTD = gammaSTD./tau;
gammaFun = cell([4,1]);
gammaFunL = zeros([4,1]);
for i = 1:4
    B = (gammaSTD(i)^2)/gammaMean(i);
    A = gammaMean(i)/B;
    y = gampdf(1:n,A,B); % mean A*B std sqrt(A*B^2)
    y = (y./max(y));
    maxL = find(y == max(y),1,'first');
    gammaFunL(i) = maxL + find(y((maxL + 1):end) < 10^(-6),1,'first');
    gammaFun{i} = y(1:gammaFunL(i));
end

%% Run the Model

% Preallocate space for the results of running the network
mV = cell([4,1]);
spkCt = cell([4,1]);
spkLoc = cell([4,1]);

% Run the two input neurons
for iNeuron = 1:2
    if strcmp(neuronTypes{iNeuron},'Ex')
        [mV{iNeuron},spkCt{iNeuron},spkLoc{iNeuron}] = simpleModel_RS(I{iNeuron},T,tau,0,memNoise(iNeuron,:));
    elseif strcmp(neuronTypes{iNeuron},'Inh')
        [mV{iNeuron},spkCt{iNeuron},spkLoc{iNeuron}] = simpleModel_FSI(I{iNeuron},T,tau,0,memNoise(iNeuron,:));
    else
        error(['Invalid neuronTypes{',num2str(iNeuron),'}'])
    end
end

% Add current pulses to input current for the other neurons based on
% the spiking activity of the input neurons
for iNeuron = 1:2
    for iSpike = find((spkLoc{iNeuron} + gammaFunL(iNeuron) - 1) <= n)
        I{3}(spkLoc{iNeuron}(iSpike):(spkLoc{iNeuron}(iSpike) + gammaFunL(iNeuron) - 1)) = ...
            I{3}(spkLoc{iNeuron}(iSpike):(spkLoc{iNeuron}(iSpike) + gammaFunL(iNeuron) - 1)) + weights(iNeuron + 2)*gammaFun{iNeuron};
        I{4}(spkLoc{iNeuron}(iSpike):(spkLoc{iNeuron}(iSpike) + gammaFunL(iNeuron) - 1)) = ...
            I{4}(spkLoc{iNeuron}(iSpike):(spkLoc{iNeuron}(iSpike) + gammaFunL(iNeuron) - 1)) + weights(iNeuron)*gammaFun{iNeuron};
    end
    for iSpike = find((spkLoc{iNeuron} + gammaFunL(iNeuron) - 1) > n)
        I{3}(spkLoc{iNeuron}(iSpike):end) = ...
            I{3}(spkLoc{iNeuron}(iSpike):end) + weights(iNeuron + 2)*gammaFun{iNeuron}(1:(n - spkLoc{iNeuron}(iSpike) + 1));
        I{4}(spkLoc{iNeuron}(iSpike):end) = ...
            I{4}(spkLoc{iNeuron}(iSpike):end) + weights(iNeuron)*gammaFun{iNeuron}(1:(n - spkLoc{iNeuron}(iSpike) + 1));
    end
end

% Run neuron 3
if strcmp(neuronTypes{3},'Ex')
    [mV{3},spkCt{3},spkLoc{3}] = simpleModel_RS(I{3},T,tau,0,memNoise(3,:));
elseif strcmp(neuronTypes{3},'Inh')
    [mV{3},spkCt{3},spkLoc{3}] = simpleModel_FSI(I{3},T,tau,0,memNoise(3,:));
else
    error('Invalid neuronTypes{3}')
end

% Add current pulses to the input current for the output neuron from neuron
% 3
iNeuron = 3;
for iSpike = find((spkLoc{iNeuron} + gammaFunL(iNeuron) - 1) <= n)
    I{4}(spkLoc{iNeuron}(iSpike):(spkLoc{iNeuron}(iSpike) + gammaFunL(iNeuron) - 1)) = ...
        I{4}(spkLoc{iNeuron}(iSpike):(spkLoc{iNeuron}(iSpike) + gammaFunL(iNeuron) - 1)) + weights(5)*gammaFun{iNeuron};
end
for iSpike = find((spkLoc{iNeuron} + gammaFunL(iNeuron) - 1) > n)
    I{4}(spkLoc{iNeuron}(iSpike):end) = ...
        I{4}(spkLoc{iNeuron}(iSpike):end) + weights(5)*gammaFun{iNeuron}(1:(n - spkLoc{iNeuron}(iSpike) + 1));
end

% Add the background current to the output neuron's input current
I{4} = I{4} + backgroundCurrent;

% Run the output neuron
if strcmp(neuronTypes{4},'Ex')
    [mV{4},spkCt{4},spkLoc{4}] = simpleModel_RS(I{4},T,tau,0,memNoise(4,:));
elseif strcmp(neuronTypes{4},'Inh')
    [mV{4},spkCt{4},spkLoc{4}] = simpleModel_FSI(I{4},T,tau,0,memNoise(4,:));
else
    error('Invalid neuronTypes{4}')
end



%% Organize the Results

% Organize the membrane potential
mVTemp = zeros([4,n]);
for iNeuron = 1:4
    mVTemp(iNeuron,:) = mV{iNeuron};
end
clear mV
mV = mVTemp;

% Organize the spike times results
spkTimes = cell([4,1]);
for iNeuron = 1:4
    spkTimes{iNeuron} = spkLoc{iNeuron}*tau;
end
clear spkLoc

% Organize the current
current = NaN([4,n]);
for iNeuron = 1:4
    current(iNeuron,:) = I{iNeuron};
end
clear I

% Organize the parameters for this model
modelParams = struct;
modelParams.noiseScale = noiseScale;
modelParams.interstimT = interstimT;
modelParams.intrastimT = intrastimT;
modelParams.nStimSets = nStimSets;
modelParams.tau = tau;
modelParams.pulseValueSTD = pulseValueSTD;
modelParams.pulseOnsetSTD = pulseOnsetSTD;
modelParams.pulseSize = pulseSize;
modelParams.weights = weights;
modelParams.gammaMean = gammaMean;
modelParams.gammaSTD = gammaSTD;
modelParams.varStim = varStim;



end