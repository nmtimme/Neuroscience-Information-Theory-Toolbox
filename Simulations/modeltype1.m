%% modeltype1 - simple stimulation of Izhikevich RS neuron
% Generates data for a simple Izhikevich regular spiking neuron. The neuron
% is stimulated with one of four types of stimuli.
%
% Syntax: [mV,spkTimes,memNoise,stim,modelParams] = modeltype1(varargin)
%
% Optional Input:
%   modeltype1(..., 'ONvsOff', pulseSize) - commands the model to run the
%       stimulus on vs. stimulus off routine. pulseSize (double) is the
%       size of the stimulus on pulse in uA. By default, the function will
%       generate data for this model with pulseSize = 500.
%   modeltype1(..., 'StimType', pulseSize) - commands the model to run a
%       routine with n stimuli. pulseSize (1 by n double array) is the
%       size of the n stimuli in uA. 
%   modeltype1(..., 'Delay', pulseSize, delay) - commands the model to run
%       the the stimulus on vs. off routine with a delay. pulseSize
%       (double) is the size of the stimulus on pulse in uA. delay (double)
%       is the delay in the stimulus reaching the neuron in milliseconds.
%   modeltype1(..., 'NonLinear', pulseSize) - commands the model to run a
%       routine with n stimuli where the stimuli are passed through a
%       non-linear filter prior to reaching the neuron. pulseSize (1 by n
%       double array) is the size of the n stimuli in uA. The nonlinear
%       filter is a quadratic function (roots = 0, 1000, max = 500).
%   modeltype1(..., 'noiseScale', noiseScale) - sets the membrane potential
%       noise strength. noiseScale (double) sets the strength (default =
%       0.1) and is the ratio of standard deviations of the noise signal. 
%       The noise is modeled as 1/f noise bound below at 10 Hz with 
%       constant power noise from 1 Hz to 10 Hz.
%   modeltype1(..., 'interstimT', interstimT) - sets the time between
%       stimulation pulses in milliseconds (double, default = 500).
%   modeltype1(..., 'intrastimT', intrastimT) - sets the duration of the
%       stimulation pulses in milliseconds (double, default = 500).
%   modeltype1(..., 'nStimSets', nStimSets) - the number of stimulation
%       sets to apply (integer, default = 50).
%   modeltype1(..., 'tau', tau) - sets the bin size or sampling in
%       milliseconds (double, default = 0.1).
%   modeltype1(..., 'pulseValueSTD', pulseValueSTD) - sets the standard
%       deviation of the Gaussian white noise added to the input current
%       received by the model neuron in pA (double, default = 5). 
%   modeltype1(..., 'pulseOnsetSTD', pulseOnsetSTD) - sets the standard
%       deviation of the Gaussian white noise added to the onset and offset
%       time of the stimuli pulses received by the model neuron in
%       milliseconds (double, default = 5).
%
% Outputs:
%   mV (double vector) - membrane potential of the model neuron in mV.
%   spkTimes (double vector) - list of times when the neuron spiked in
%       milliseconds.
%   memNoise (double vector) - the membrane noise of the model neuron in
%       mV.
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
%       [mV,spkTimes,memNoise,stim,modelParams] = modeltype1;
%   Three Types of Stimulation Strengths
%       [mV,spkTimes,memNoise,stim,modelParams] = modeltype1('StimType', [200,500,800]);
%   Weaker Stimulation with a 100 ms Delay and No Stimulus Current Noise
%       [mV,spkTimes,memNoise,stim,modelParams] = modeltype1('Delay', 500, 100, 'pulseValueSTD', 0, 'pulseOnsetSTD', 0);
%   
% Other m-files required: inverseFNoise, simpleModel_RS
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% February 2016; Last revision: 8-Feb-2016


function [mV,spkTimes,memNoise,stim,modelParams] = modeltype1(varargin)
%% Parse command line for parameters


% Set default values
noiseScale = 0.1;
interstimT = 500;
intrastimT = 500;
nStimSets = 50;
tau = 0.1;
pulseValueSTD = 5;
pulseOnsetSTD = 5;
stimFormat = 'ONvsOFF';
pulseSize = 500;
delay = NaN;



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
        case 'ONvsOFF',             stimFormat = varargin{iVarArg}; pulseSize = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'StimType',            stimFormat = varargin{iVarArg}; pulseSize = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'Delay',               stimFormat = varargin{iVarArg}; pulseSize = varargin{iVarArg+1}; delay = varargin{iVarArg+2}; iVarArg = iVarArg + 2;
        case 'NonLinear',           stimFormat = varargin{iVarArg}; pulseSize = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(MODELTYPE1) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Do Some Preprocessing

% Add zero stimulation where appropriate
if strcmp(stimFormat,'ONvsOFF') || strcmp(stimFormat,'Delay')
    pulseSize = [pulseSize,0];
end
nPulses = length(pulseSize);

% Find the total length of the recording in milliseconds
T = nStimSets*nPulses*intrastimT + nStimSets*nPulses*interstimT + interstimT;

% Find the total number of time bins
n = round(T/tau);

%% Make the Input Current, Add Noise, and Make the Membrane Noise

% Make the input current
vInput = [];
stims = [];
for iPulse = 1:nPulses
    vInput = [vInput,zeros([1,round(interstimT/tau)]),pulseSize(iPulse)*ones([1,round(intrastimT/tau)])];
    stims = [stims,zeros([1,round(interstimT/tau)]),iPulse,zeros([1,round(intrastimT/tau) - 1])];
end
vInput = repmat(vInput,[1,nStimSets]);
vInput = [vInput,zeros([1,round(interstimT/tau)])];
while length(vInput) > n
    vInput(end) = [];
end
while length(vInput) < n
    vInput = [vInput,0];
end
stims = repmat(stims,[1,nStimSets]);
stim = struct;
stim.times = tau*find(stims ~= 0);
stim.types = stims(stims ~= 0);

% Keep a copy of the original stimulation current
stim.currentOriginal = vInput;

% Add noise to the pulse onset and offset
transpoints = find(diff(vInput) ~= 0);
noise = round((pulseOnsetSTD/tau)*randn(size(transpoints)));
for itrans = 1:length(transpoints)
    if noise(itrans) > 0
        vInput((transpoints(itrans) + 1):(transpoints(itrans) + noise(itrans))) = vInput(transpoints(itrans));
    elseif noise(itrans) < 0
        vInput((transpoints(itrans) + noise(itrans) + 1):transpoints(itrans)) = vInput(transpoints(itrans) + 1);
    end
end

% Add noise to the input current
vInput = vInput + pulseValueSTD*randn(size(vInput));

% If this is the non-linear stimulation model, transform the stimulation
if strcmp(stimFormat,'NonLinear')
    % Quadratic Polynomial (roots = 0, 1000) (max = 500)
    vInput = -(1/500)*vInput.*(vInput - 1000);
end

% If this is the delay stimulation model, add a delay to the input current
% the neuron will see
if strcmp(stimFormat,'Delay')
    vInput = [pulseValueSTD*randn([1,round(delay/tau)]),vInput];
    vInput((end - round(delay/tau) + 1):end) = [];
end

% Keep a copy of the stimulation current received by the neuron
stim.currentRec = vInput;

% Make the membrane noise
nzInput = inverseFNoise(n, 1, 1000/tau, 'minF', 10, 'constLow', 1);
nzInput = nzInput./std(nzInput);
nzInput = nzInput*noiseScale;

% Keep a copy of the membrane noise
memNoise = nzInput;

%% Run the Model and Organize the Results

% Run the model
[mV,spkCt_E1,spkTimes]=simpleModel_RS(vInput,T,tau,0,nzInput);

% Convert the spike location to time in milliseconds
spkTimes = spkTimes*tau;

% Organize the parameters for this model
modelParams = struct;
modelParams.noiseScale = noiseScale;
modelParams.interstimT = interstimT;
modelParams.intrastimT = intrastimT;
modelParams.nStimSets = nStimSets;
modelParams.tau = tau;
modelParams.pulseValueSTD = pulseValueSTD;
modelParams.pulseOnsetSTD = pulseOnsetSTD;
modelParams.stimFormat = stimFormat;
modelParams.pulseSize = pulseSize;
modelParams.delay = delay;

end