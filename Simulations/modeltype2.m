%% modeltype2 - three Izhikevich neuron system (2 RS, 1 FS)
% Generates data for a three Izhikevich neuron network. Pulse stimulus is
% applied to a RS neuron (E1) that sends connections to another RS neuron 
% (E2) and an inhibitory FS neuron (I1). I1 also sends a connection to E2.
% The stimulus to E1 can either be a uniform square wave (with noise) or
% sequences of instantaneous current pulses.
%
% Syntax: [mV,spkTimes,memNoise,current,stim,modelParams] = modeltype2(pulseSize, weights, varargin)
%
% Inputs:
%   pulseSize - sets the sizes in pA of the n stimulus pulses for square
%       waves or the spike rate in Hz for the spike pulses (1 by n double
%       vector). (See stimType option below.)
%   weights - sets the maximum current in pA supplied from presynaptic
%       spikes to postsynaptic neurons (1 by 3 double vector). weights(1)
%       is the weight from E1 to E2. weights(2) is the weight from E1 to
%       I1. weights(3) is the weight from I1 to E2 (it should be
%       negative). The weight value corresponds to the maximum current
%       supplied using a gamma pdf.
% 
% Optional Inputs:
%   modeltype2(..., 'stimType', stimType) - sets the type of stimulation to
%       apply (string). Options are:
%         stimType = 'UniSqPulse': a uniform square pulse is applied
%           (default). pulseSize controls the size of the pulse in pA. 
%         stimType = 'UnqSpkPulse': instantaneous spike pulses of magnitude
%           spkCurrent are applied. pulseSize control the spike rate in Hz
%           of the spike pulses. Each sequence is unique.
%         stimType = 'RepSpkPulse': instantaneous spike pulses of magnitude
%           spkCurrent are applied. pulseSize control the spike rate in Hz
%           of the spike pulses. The spike pulses are repeated across
%           stimulation sets.
%   modeltype2(..., 'noiseScale', noiseScale) - sets the membrane potential
%       noise strength. noiseScale (1 by 3 double vector) sets the strength 
%       (default = 0.1 for all neurons) and is the ratio of standard 
%       deviations of the noise signal.  The noise is modeled as 1/f noise 
%       bound below at 10 Hz with constant power noise from 1 Hz to 10 Hz. 
%   modeltype2(..., 'interstimT', interstimT) - sets the time between
%       stimulation pulses in milliseconds (double, default = 500).
%   modeltype2(..., 'intrastimT', intrastimT) - sets the duration of the
%       stimulation pulses in milliseconds (double, default = 500). If
%       spike stimuli are applied, this value sets the window over which
%       the spike stimuli can occur.
%   modeltype2(..., 'nStimSets', nStimSets) - the number of stimulation
%       sets to apply (integer, default = 50).
%   modeltype2(..., 'tau', tau) - sets the bin size or sampling in
%       milliseconds (double, default = 0.1).
%   modeltype2(..., 'pulseValueSTD', pulseValueSTD) - sets the standard
%       deviation of the Gaussian white noise added to the input current
%       received by the model neuron in pA (double, default = 5). 
%   modeltype2(..., 'pulseOnsetSTD', pulseOnsetSTD) - sets the standard
%       deviation of the Gaussian white noise added to the onset and offset
%       time of the stimuli pulses received by the model neuron in
%       milliseconds (double, default = 5). This option only applies to
%       square pulses (stimType = 'UniSqPulse').
%   modeltype2(..., 'gammaMean', gammaMean) - sets the mean of the gamma
%       pdf used to model current injection to postsynaptic neurons from
%       presynaptic spikes (1 by 3 double vector, default = [30,30,30]). 
%       gammaMean(1) applies to the connection from E1 to E2. gammaMean(2) 
%       applies to the connection from E1 to I1. gammaMean(3) applies to
%       the connection from I1 to E2. 
%   modeltype2(..., 'gammaSTD', gammaSTD) - sets the standard deviation of
%       the gamma pdf used to model current injection to postsynaptic
%       neurons from presynaptic spikes (1 by 3 double vector, default =
%       [20,20,20]). gammaSTD(1) applies to the connection from E1 to E2.
%       gammaSTD(2) applies to the connection from E1 to I1. gammaSTD(3)
%       applies to the connection from I1 to E2.
%   modeltype2(..., 'spkCurrent', spkCurrent) - sets the magnitude of the
%       spike current pulses in pA (double, default = 500). This parameter
%       does not apply to square pulses.
%   modeltype2(..., 'spkDur', spkDur) - sets the duration of the spike
%       current stimulus pulses in milliseconds(double, default = 5). This
%       parameter does not apply to square pulses.
%
% Outputs:
%   mV (3 by number of time steps double array) - membrane potential of the 
%       model neurons in mV. Row 1 is E1, row 2 is I1, row 3 is E2.
%   spkTimes (3 by 1 cell array, each element is a double vector) - list of 
%       times when a neuron spiked in milliseconds. spkTimes{1} is E1,
%       spkTimes{2} is I1, spkTimes{3} is E2.
%   memNoise (3 by number of time steps double vector) - the membrane noise 
%       of the model neurons in mV. Row 1 is E1, row 2 is I1, row 3 is E2.
%   current (3 by number of time steps double vector) - the input current
%       to each neuron in pA. Row 1 is E1 (and is identical to
%       stim.currentRec). Row 2 is I1. Row 3 is E2.
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
%   Stimulation On Vs. Off Data
%       [mV,spkTimes,memNoise,stim,modelParams] = modeltype2([500,0],[100,100,100]);
%   No Inhibitory Interaction
%       [mV,spkTimes,memNoise,stim,modelParams] = modeltype2([500,0],[100,100,0]);
%   Three Different Stimulation Levels
%       [mV,spkTimes,memNoise,stim,modelParams] = modeltype2([500,250,0],[100,100,100]);
%   
% Other m-files required: inverseFNoise, simpleModel_RS, simpleModel_FSI
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% February 2016; Last revision: 19-Feb-2016


function [mV,spkTimes,memNoise,current,stim,modelParams] = modeltype2(pulseSize, weights, varargin)
%% Parse command line for parameters


% Set default values
noiseScale = [0.1,0.1,0.1];
interstimT = 500;
intrastimT = 500;
nStimSets = 50;
tau = 0.1;
pulseValueSTD = 5;
pulseOnsetSTD = 5;
gammaMean = [30,30,30];
gammaSTD = [20,20,20];
stimType = 'UniSqPulse';
spkCurrent = 500;
spkDur = 5;



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
        case 'stimType',            stimType = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'spkCurrent',          spkCurrent = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'spkDur',              spkDur = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(MODELTYPE2) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Do Some Preprocessing

% Error check the inhibitory weight
if weights(3) > 0
    error('The inhibitory weight (weights(3)) should be less than or equal to zero.')
end

% Error check the pulse spike duration
if round(spkDur/tau) < 1
    error('Stimulus spike pulses are too short.')
end

% Find the number of pulses in the stimulation
nPulses = length(pulseSize);

% Find the total length of the recording in milliseconds
T = nStimSets*nPulses*intrastimT + nStimSets*nPulses*interstimT + interstimT;

% Find the total number of time bins
n = round(T/tau);

%% Make the Input Current, Add Noise, and Make the Membrane Noise

% Make the input current to E1
if strcmp(stimType,'UniSqPulse')
    
    % Make uniform square pulses
    
    vInputE1 = [];
    stims = [];
    for iPulse = 1:nPulses
        vInputE1 = [vInputE1,zeros([1,round(interstimT/tau)]),pulseSize(iPulse)*ones([1,round(intrastimT/tau)])];
        stims = [stims,zeros([1,round(interstimT/tau)]),iPulse,zeros([1,round(intrastimT/tau) - 1])];
    end
    vInputE1 = repmat(vInputE1,[1,nStimSets]);
    vInputE1 = [vInputE1,zeros([1,round(interstimT/tau)])];
    while length(vInputE1) > n
        vInputE1(end) = [];
    end
    while length(vInputE1) < n
        vInputE1 = [vInputE1,0];
    end
    stims = repmat(stims,[1,nStimSets]);
    stim = struct;
    stim.times = tau*find(stims ~= 0);
    stim.types = stims(stims ~= 0);
    
    % Keep a copy of the original stimulation current
    stim.currentOriginal = vInputE1;
    
    % Add noise to the pulse onset and offset
    transpoints = find(diff(vInputE1) ~= 0);
    noise = round((pulseOnsetSTD/tau)*randn(size(transpoints)));
    for itrans = 1:length(transpoints)
        if noise(itrans) > 0
            vInputE1((transpoints(itrans) + 1):(transpoints(itrans) + noise(itrans))) = vInputE1(transpoints(itrans));
        elseif noise(itrans) < 0
            vInputE1((transpoints(itrans) + noise(itrans) + 1):transpoints(itrans)) = vInputE1(transpoints(itrans) + 1);
        end
    end
    
elseif strcmp(stimType,'UnqSpkPulse')
    
    % Make unique spike pulses
    
    % Convert pulse rates to probabilities
    pulseSize = pulseSize./(1000/tau);
    
    vInputE1 = [];
    stims = [];
    for iPulse = 1:nPulses
        vInputE1 = [vInputE1,zeros([1,round(interstimT/tau)]),pulseSize(iPulse)*ones([1,round(intrastimT/tau)])];
        stims = [stims,zeros([1,round(interstimT/tau)]),iPulse,zeros([1,round(intrastimT/tau) - 1])];
    end
    vInputE1 = repmat(vInputE1,[1,nStimSets]);
    vInputE1 = [vInputE1,zeros([1,round(interstimT/tau)])];
    while length(vInputE1) > n
        vInputE1(end) = [];
    end
    while length(vInputE1) < n
        vInputE1 = [vInputE1,0];
    end
    stims = repmat(stims,[1,nStimSets]);
    stim = struct;
    stim.times = tau*find(stims ~= 0);
    stim.types = stims(stims ~= 0);
    
    % Convert from probabilities to pulse spikes
    vInputE1(rand(size(vInputE1)) < vInputE1) = 1;
    vInputE1(vInputE1 < 1) = 0;
    
    % Include the pulse spike durations
    spkPulseLoc = find(vInputE1 == 1);
    for ispk = 1:length(spkPulseLoc)
        vInputE1(spkPulseLoc(ispk):(spkPulseLoc(ispk) + round(spkDur/tau) - 1)) = spkCurrent;
    end
    
    % Keep a copy of the original stimulation current
    stim.currentOriginal = vInputE1;
    
    
elseif strcmp(stimType,'RepSpkPulse')
    
    % Make repeated spike pulses
    
    % Convert pulse rates to probabilities
    pulseSize = pulseSize./(1000/tau);
    
    vInputE1 = [];
    stims = [];
    for iPulse = 1:nPulses
        vInputE1 = [vInputE1,zeros([1,round(interstimT/tau)]),pulseSize(iPulse)*ones([1,round(intrastimT/tau)])];
        stims = [stims,zeros([1,round(interstimT/tau)]),iPulse,zeros([1,round(intrastimT/tau) - 1])];
    end
    
    % Convert from probabilities to pulse spikes
    vInputE1(rand(size(vInputE1)) < vInputE1) = 1;
    vInputE1(vInputE1 < 1) = 0;
    
    % Include the pulse spike durations
    spkPulseLoc = find(vInputE1 == 1);
    for ispk = 1:length(spkPulseLoc)
        vInputE1(spkPulseLoc(ispk):(spkPulseLoc(ispk) + round(spkDur/tau) - 1)) = spkCurrent;
    end
    
    % Duplicate the spike sequence
    vInputE1 = repmat(vInputE1,[1,nStimSets]);
    vInputE1 = [vInputE1,zeros([1,round(interstimT/tau)])];
    while length(vInputE1) > n
        vInputE1(end) = [];
    end
    while length(vInputE1) < n
        vInputE1 = [vInputE1,0];
    end
    stims = repmat(stims,[1,nStimSets]);
    stim = struct;
    stim.times = tau*find(stims ~= 0);
    stim.types = stims(stims ~= 0);
    
    % Keep a copy of the original stimulation current
    stim.currentOriginal = vInputE1;
    
end

% Add noise to the input current
vInputE1 = vInputE1 + pulseValueSTD*randn(size(vInputE1));

% Keep a copy of the stimulation current received by the neuron
stim.currentRec = vInputE1;

% Make the membrane noise
memNoise = NaN([3,n]);
for i = 1:3
    memNoise(i,:) = inverseFNoise(n, 1, 1000/tau, 'minF', 10, 'constLow', 1);
    memNoise(i,:) = memNoise(i,:)./std(memNoise(i,:));
    memNoise(i,:) = memNoise(i,:)*noiseScale(i);
end

% Make the necessary gamma pdfs for the current injected by spikes
gammaMean = gammaMean./tau;
gammaSTD = gammaSTD./tau;
gammaFun = cell([3,1]);
gammaFunL = zeros([3,1]);
for i = 1:3
    B = (gammaSTD(i)^2)/gammaMean(i);
    A = gammaMean(i)/B;
    y = gampdf(1:n,A,B); % mean A*B std sqrt(A*B^2)
    y = (y./max(y));
    maxL = find(y == max(y),1,'first');
    gammaFunL(i) = maxL + find(y((maxL + 1):end) < 10^(-6),1,'first');
    gammaFun{i} = weights(i)*y(1:gammaFunL(i)); % Include weight here
end

%% Run the Model

% Run E1
[mV_E1,spkCt_E1,spkLoc_E1] = simpleModel_RS(vInputE1,T,tau,0,memNoise(1,:));

% Make the current received by E2 and I1 from E1
if spkCt_E1 >= 1
    vInputE1E2 = zeros(1,n);
    vInputE1I1 = zeros(1,n);
    for i=find((spkLoc_E1 + gammaFunL(1) - 1) <= n)
        vInputE1E2(spkLoc_E1(i):(spkLoc_E1(i) + gammaFunL(1) - 1)) = vInputE1E2(spkLoc_E1(i):(spkLoc_E1(i) + gammaFunL(1) - 1)) + gammaFun{1};
    end
    for i=find((spkLoc_E1 + gammaFunL(1) - 1) > n)
        vInputE1E2(spkLoc_E1(i):end) = vInputE1E2(spkLoc_E1(i):end) + gammaFun{1}(1:(n - spkLoc_E1(i) + 1));
    end
    for i=find((spkLoc_E1 + gammaFunL(2) - 1) <= n)
        vInputE1I1(spkLoc_E1(i):(spkLoc_E1(i) + gammaFunL(2) - 1)) = vInputE1I1(spkLoc_E1(i):(spkLoc_E1(i) + gammaFunL(2) - 1)) + gammaFun{2};
    end
    for i=find((spkLoc_E1 + gammaFunL(2) - 1) > n)
        vInputE1I1(spkLoc_E1(i):end) = vInputE1I1(spkLoc_E1(i):end) + gammaFun{2}(1:(n - spkLoc_E1(i) + 1));
    end
else
    vInputE1E2 = zeros([1,n]);
    vInputE1I1 = zeros([1,n]);
end


% Run I1
[mV_I1,spkCt_I1,spkLoc_I1] = simpleModel_FSI(vInputE1I1,T,tau,0,memNoise(2,:));

% Make the current received by E2 from I1
if spkCt_I1 >= 1
    vInputI1E2 = zeros(1,n);
    for i=find((spkLoc_I1 + gammaFunL(3) - 1) <= n)
        vInputI1E2(spkLoc_I1(i):(spkLoc_I1(i) + gammaFunL(3) - 1)) = vInputI1E2(spkLoc_I1(i):(spkLoc_I1(i) + gammaFunL(3) - 1)) + gammaFun{3};
    end
    for i=find((spkLoc_I1 + gammaFunL(3) - 1) > n)
        vInputI1E2(spkLoc_I1(i):end) = vInputI1E2(spkLoc_I1(i):end) + gammaFun{3}(1:(n - spkLoc_I1(i) + 1));
    end
else
    vInputI1E2 = zeros([1,n]);
end

% Sum the currents
vInputE2 = vInputI1E2 + vInputE1E2;

% Run E2
[mV_E2,spkCt_E2,spkLoc_E2] = simpleModel_RS(vInputE2,T,tau,0,memNoise(3,:));




%% Organize the Results

% Organize the membrane potential results
mV = NaN([3,n]);
mV(1,:) = mV_E1;
mV(2,:) = mV_I1;
mV(3,:) = mV_E2;

% Organize the spike times results
spkTimes = cell([3,1]);
spkTimes{1} = spkLoc_E1*tau;
spkTimes{2} = spkLoc_I1*tau;
spkTimes{3} = spkLoc_E2*tau;

% Organize the current
current = NaN([3,n]);
current(1,:) = vInputE1;
current(2,:) = vInputE1I1;
current(3,:) = vInputE2;

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



end