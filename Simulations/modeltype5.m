%% modeltype5 - model canonical experiments
%
% Syntax: [data] = modeltype5(expType, varargin)
%
% Inputs:
%   expType - sets the type of experiment to run (string). Options are:
%           'CenterOut': a movement direction encoding by primary motor 
%               cortex neurons experiment similar to Georgopoulis. 
%           'CenterSurround': visual stimulus encoding by center surround
%               retinal ganglion cell
%           'PlaceCells': neurons fire preferentially based on the animal's
%               location in a 2-D space similarly to place cells
%           'SensoryHabit': synapse between a model sensory neuron and a
%               model motor neuron weakens with use (habituation), similar 
%               to simple models of gill-withdrawl reflex in Aplysia.
%
% Citations
%
% - Georgopoulos, A.P., Kalaska, J.F., Caminiti, R., Massey, J.T., 1982. On
%   the relations between the direction of two-dimensional arm movements 
%   and cell discharge in primate motor cortex. Journal of Neuroscience 2, 
%   1527-1537.
% - Castellucci, V., Pinsker, H., Kupfermann, I., Kandel, E.R., 1970.
%   Neuronal mechanisms of habituation and dishabituation of the 
%   gill-withdrawl reflex in aplysia. Science 167, 1745-1748.
% - Kuffler, S.W., 1953. Discharge patterns and functional organization of
%   mammalian retina. Journal of Neurophysiology 16, 37-68.
% - O'Keefe, J., Dostrovsky, J., 1971. The hippocampus as a spatial map.
%   Brain Research 34, 171-175.
% - Bear, M.F., Connors, B.W., Paradiso, M.A., 2007. Neuroscience: 
%   exploring the brain, Third ed. Lippincott Williams and Wilkins, 
%   Baltimore, Maryland.
% 
% Optional Inputs:
%   modeltype5(..., 'parameters', parameters) - sets experiment specific
%       parameters (structure). The precise format of parameters is
%       dependent upon expType. Options are:
%           expType = 'CenterOut': 
%           expType = 'CenterSurround':
%           expType = 'PlaceCells':
%           expType = 'SensoryHabit':
%
% Outputs:
%   data (structure): contains the results of the model. The precise format
%       of data is dependent upon expType. Options are:
%           expType = 'CenterOut': 
%           expType = 'CenterSurround':
%           expType = 'PlaceCells':
%           expType = 'SensoryHabit':
%
% Examples:
%   
% Other m-files required: inverseFNoise, simpleModel_RS, simpleModel_FSI
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% March 2016; Last revision: 24-Mar-2016


function [data] = modeltype5(expType, varargin)
%% Parse command line for parameters


% Set default parameters

% Center Out Model
parametersCO = struct;
parametersCO.maxFR = 100;
parametersCO.meanLead = -100;
parametersCO.stdLead = 10;
parametersCO.stdMP = 100;
parametersCO.nAngles = 8;
parametersCO.nNeurons = 10;
parametersCO.nTrials = 150;
parametersCO.intertrialT = 1000;
parametersCO.respRange = [0,1];

% Center Surround parameters
parametersCS = struct;
parametersCS.smallR = 0.08;
parametersCS.largeR = 0.16;
parametersCS.bgFR = 30;
parametersCS.centerFR = 100;
parametersCS.surroundFR = 1;
parametersCS.nStims = 300;
parametersCS.interstimT = 100;
parametersCS.intrastimT = 500;
parametersCS.nNeurons = 200;

% Place Cells parameters
parametersPC = struct;
parametersPC.nNeurons = 10;
parametersPC.maxFR = 100;
parametersPC.bgFR = 20;
parametersPC.spatialSTD = 0.2;
parametersPC.stepDur = 10;
parametersPC.nSteps = 20000;
parametersPC.nStepsDim = 20;
parametersPC.respRange = [0,1];

% Sensory Habituation parameters
parametersSH = struct;
parametersSH.nStimSets = 50;
parametersSH.interstimT = 500;
parametersSH.intrastimT = 500;
parametersSH.pulseSize = 500;
parametersSH.maxI = 200;
parametersSH.minI = 30;
parametersSH.tau = 0.1;
parametersSH.pulseOnsetSTD = 5;
parametersSH.pulseValueSTD = 5;
parametersSH.gammaMean = 30*ones([1,4]);
parametersSH.gammaSTD = 20*ones([1,4]);
parametersSH.noiseScale = [0.1,0.1];




iVarArg = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'parameters',          parametersTemp = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(MODELTYPE5) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Do Some Preprocessing

% Identify the parameters and error check
if strcmp(expType,'CenterOut')
    
    try 
        parameters = parametersTemp;
    catch
        parameters = parametersCO;
    end
    
    % Make sure there are enough angles
    if parameters.nAngles < 2
        error('Not enough angles')
    end
    
    % Make sure the max firing rate can be capture by this bin size
    if parameters.maxFR > 1000
        error('The maximum firing rate is too high for this bin size (1 ms).')
    end
    
    % Make sure time between trials is at least 5 standard deviations of
    % the modification profile to make sure trials do not overlap
    if parameters.stdMP > (parameters.intertrialT/5)
        error('Inter trial time is too small for this modification profile.')
    end
    
    
elseif strcmp(expType,'CenterSurround')
    
    try 
        parameters = parametersTemp;
    catch
        parameters = parametersCS;
    end
    
    % Make sure the max firing rate can be capture by this bin size
    if parameters.centerFR > 1000
        error('The center region firing rate is too high for this bin size (1 ms).')
    end
    
    
elseif strcmp(expType,'PlaceCells')
    
    try 
        parameters = parametersTemp;
    catch
        parameters = parametersPC;
    end
    
    
elseif strcmp(expType,'SensoryHabit')
    
    try 
        parameters = parametersTemp;
    catch
        parameters = parametersSH;
    end
    
end

% Preallocate the data structure
data = struct;

%% Run the Model

if strcmp(expType,'CenterOut') 
    
    %% Center Out Model
    
    % Randomly select preferential angles
    prefA = 2*pi*rand([1,parameters.nNeurons]);
    
    % Randomly select the strength of the response
    respS = (parameters.respRange(2) - parameters.respRange(1))*rand([1,parameters.nNeurons]) + parameters.respRange(1);
    
    % Figure out the background spike probability per bin based on the
    % background firing rate
    pBG = (parameters.maxFR/2) / 1000;
    
    % Randomly select the directions of movement
    angles = linspace(0,2*pi,parameters.nAngles + 1);
    angles(end) = [];
    goalA = angles(randi(parameters.nAngles,[1,parameters.nTrials]));
    
    % Figure out how many total bins will be in the recording
    nTime = parameters.intertrialT*(parameters.nTrials + 1);
    
    % Make a list of the trial times
    trialTimes = (parameters.intertrialT + 1):parameters.intertrialT:nTime;
    
    % Loop through each neuron
    data.spkTimes = cell([parameters.nNeurons,1]);
    for iNeuron = 1:parameters.nNeurons
        
        % Make a vector of the spiking probability for this bin
        p = pBG*ones([1,nTime]);
        
        % Randomly select the lead time
        leadTime = parameters.stdLead*randn([1,parameters.nTrials]) + parameters.meanLead;
        
        for iTrial = 1:parameters.nTrials
            
            % Make a list of the time values relative to this trial time
            time = (1:nTime) - trialTimes(iTrial);
            
            % Modify the firing probability based on this trial
            p = p + pBG*respS(iNeuron)*cos(goalA(iTrial) - prefA(iNeuron))*exp(-((time - leadTime(iTrial)).^2)/(2*parameters.stdMP^2));
            
        end
        
        % Randomly insert firing
        data.spkTimes{iNeuron} = find(p > rand([1,nTime]));
        
    end
    
    % Record other relevant values
    data.prefA = prefA;
    data.respS = respS;
    data.goalA = goalA;
    data.angles = angles;
    data.trialTimes = trialTimes;
    data.parameters = parameters;
    data.expType = expType;
            
            
        
elseif strcmp(expType,'CenterSurround')
    
    %% Center Surround Model
    
    % Randomly select the cell locations, except for the first neuron,
    % which goes at the center
    cellLocs = [0.5,0.5;rand([parameters.nNeurons - 1,2])];
    
    % Make the mirror cell preference locations for periodic boundary
    % conditions
    temp = cell([8,1]);
    iMirror = 1;
    for i = [-1,0,1]
        for j = [-1,0,1]
            if (i ~= 0) || (j ~= 0)
                temp{iMirror} = cellLocs + repmat([i,j],[parameters.nNeurons,1]);
                iMirror = iMirror + 1;
            end
        end
    end            
    mirrorCellLocs = cat(3,temp{1},temp{2},temp{3},temp{4},temp{5},temp{6},temp{7},temp{8});
    
    % Randomly select the stimulus locations
    stimLocs = rand([parameters.nStims,2]);
    %stimLocs = repmat(rand([1,2]),[parameters.nStims,1]); % For testing
    
    % Figure out how many total bins will be in the recording
    nTime = parameters.interstimT*(parameters.nStims + 1) + parameters.intrastimT*parameters.nStims;
    
    % Find the times for the stimulus onset
    stimTimes = (parameters.interstimT + 1):(parameters.interstimT + parameters.intrastimT):nTime;
    
    % Figure out the spiking probabilities based on the firing rates
    pBG = (parameters.bgFR) / 1000;
    pCenter = (parameters.centerFR) / 1000;
    pSurround = (parameters.surroundFR) / 1000;
    
    % Loop through each neuron
    data.spkTimes = cell([parameters.nNeurons,1]);
    for iNeuron = 1:parameters.nNeurons
        
        % Make a vector of the spiking probability for this bin
        p = pBG*ones([1,nTime]);
        
        for iStim = 1:parameters.nStims
            
            % Calculate the minimum distance to the center of the receptive
            % field
            r = min(sqrt(sum(([cellLocs(iNeuron,:);squeeze(permute(mirrorCellLocs(iNeuron,:,:),[3,2,1]))] - repmat(stimLocs(iStim,:),[9,1])).^2,2)));
            
            % Check if this stimulus is in the surround or the center
            if r <= parameters.smallR
                p(stimTimes(iStim):(stimTimes(iStim) + parameters.intrastimT - 1)) = pCenter;
            elseif r <= parameters.largeR
                p(stimTimes(iStim):(stimTimes(iStim) + parameters.intrastimT - 1)) = pSurround;
            end
            
        end
        
        % Randomly insert firing
        data.spkTimes{iNeuron} = find(p > rand([1,nTime]));
        
    end
    
    % Record other relevant values
    data.cellLocs = cellLocs;
    data.stimLocs = stimLocs;
    data.stimTimes = stimTimes;
    data.parameters = parameters;
    data.expType = expType;
    
    
    
elseif strcmp(expType,'PlaceCells')
    
    %% Place Cells Model
    
    % Randomly select the cell preferential firing locations, except for
    % the first neuron, which we force to be in the middle
    cellPrefLocs = [0.5,0.5;rand([parameters.nNeurons - 1,2])];
    
    % Make the mirror cell preference locations for periodic boundary
    % conditions
    temp = cell([8,1]);
    iMirror = 1;
    for i = [-1,0,1]
        for j = [-1,0,1]
            if (i ~= 0) || (j ~= 0)
                temp{iMirror} = cellPrefLocs + repmat([i,j],[parameters.nNeurons,1]);
                iMirror = iMirror + 1;
            end
        end
    end            
    mirrorCellPrefLocs = cat(3,temp{1},temp{2},temp{3},temp{4},temp{5},temp{6},temp{7},temp{8});
    
    % Randomly select the strength of the response
    respS = (parameters.respRange(2) - parameters.respRange(1))*rand([1,parameters.nNeurons]) + parameters.respRange(1);
    
    % Figure out how many total bins will be in the recording
    nTime = parameters.stepDur*parameters.nSteps;
    
    % Figure out the allowed animal locations
    allowedLocs = linspace(0,1,parameters.nStepsDim);
    
    % Perform the random walk
    animalLoc = zeros([parameters.nSteps,2]);
    animalLoc(1,:) = round(parameters.nStepsDim/2);
    temp1 = [parameters.nStepsDim,1,2];
    temp2 = [parameters.nStepsDim - 1,parameters.nStepsDim,1];
    temp3 = [-1,0,1];
    for iStep = 2:parameters.nSteps
        for iDim = 1:2
            if animalLoc(iStep - 1,iDim) == 1
                animalLoc(iStep,iDim) = temp1(randi(3));
            elseif animalLoc(iStep - 1,iDim) == parameters.nStepsDim
                animalLoc(iStep,iDim) = temp2(randi(3));
            else
                animalLoc(iStep,iDim) = animalLoc(iStep - 1,iDim) + temp3(randi(3));
            end
        end
    end
    
    % Expand the random walk to note that each step takes more than 1 bin
    tempx = repmat(animalLoc(:,1)',[parameters.stepDur,1]);
    tempy = repmat(animalLoc(:,2)',[parameters.stepDur,1]);
    animalLoc = [tempx(:),tempy(:)];
    
    % Convert steps to physical location
    animalLoc = allowedLocs(animalLoc);
    
    % Figure out the spiking probabilities based on the firing rates
    pBG = (parameters.bgFR) / 1000;
    pMax = (parameters.maxFR) / 1000;
    
    % Loop through each neuron
    data.spkTimes = cell([parameters.nNeurons,1]);
    for iNeuron = 1:parameters.nNeurons
        
        % Make a vector of the spiking probability for this bin
        p = pBG*ones([nTime,1]);
        
        % Calculate the distance between the animal and the preferred
        % firing location for all time
        r = min(sqrt(sum((cat(3,repmat(cellPrefLocs(iNeuron,:),[nTime,1]),repmat(mirrorCellPrefLocs(iNeuron,:,:),[nTime,1,1])) - repmat(animalLoc,[1,1,9])).^2,2)),[],3);
        
        % Calculate the increase in firing probability based on the
        % distance to the preferred location
        p = p + respS(iNeuron)*(pMax - pBG)*exp(-(r.^2)/(2*parameters.spatialSTD^2));
        
        % Randomly insert firing
        data.spkTimes{iNeuron} = find(p > rand([nTime,1]));
        
    end
    
    % Record other relevant values
    data.cellPrefLocs = cellPrefLocs;
    data.animalLoc = animalLoc;
    data.respS = respS;
    data.parameters = parameters;
    data.expType = expType;
    
    
    
    
elseif strcmp(expType,'SensoryHabit')
    
    %% Sensory Habituation Model
    
    % Make the stimulus, add noise, and make the membrane noise
    
    % Make the pulse size vector
    pulseSize = [0,parameters.pulseSize];
    
    % Find the number of pulses in one stimulation set
    nPulses = size(pulseSize,2);
    
    % Find the total length of the recording in milliseconds
    T = parameters.nStimSets*nPulses*parameters.intrastimT + parameters.nStimSets*nPulses*parameters.interstimT + parameters.interstimT;
    
    % Find the total number of time bins
    n = round(T/parameters.tau);
    
    % Figure out the stimulus times
    stimTimes = (parameters.interstimT + parameters.tau):(parameters.interstimT + parameters.intrastimT):T;
    nStims = length(stimTimes);
    
    % Make the stimulus currents
    I = cell([2,1]);
    for iNeuron = 1:2
        I{iNeuron} = zeros([1,n]);
    end
    time = (1:n)*parameters.tau;
    stimTypes = zeros([1,nStims]);
    for iStim = 1:nStims
        if rem(iStim,2) == 1
            I{1}((time >= stimTimes(iStim)) & (time < (stimTimes(iStim) + parameters.intrastimT))) = pulseSize(1);
            stimTypes(iStim) = 1;
        else
            I{1}((time >= stimTimes(iStim)) & (time < (stimTimes(iStim) + parameters.intrastimT))) = pulseSize(2);
            stimTypes(iStim) = 2;
        end
    end
    
    % Keep a copy of the original stimulation current
    currentOriginal = I{1};
    
    % Add noise to the pulse onset and offset
    transpoints = find(diff(I{1}) ~= 0);
    noise = round((parameters.pulseOnsetSTD/parameters.tau)*randn(size(transpoints)));
    for itrans = 1:length(transpoints)
        if noise(itrans) > 0
            I{1}((transpoints(itrans) + 1):(transpoints(itrans) + noise(itrans))) = I{1}(transpoints(itrans));
        elseif noise(itrans) < 0
            I{1}((transpoints(itrans) + noise(itrans) + 1):transpoints(itrans)) = I{1}(transpoints(itrans) + 1);
        end
    end
    
    % Add noise to the input current
    I{1} = I{1} + parameters.pulseValueSTD*randn(size(I{1}));
    
    % Keep a copy of the stimulation current received by the neuron
    currentRec = I{1};
    
    % Make the membrane noise
    memNoise = NaN([2,n]);
    for i = 1:2
        memNoise(i,:) = inverseFNoise(n, 1, 1000/parameters.tau, 'minF', 10, 'constLow', 1);
        memNoise(i,:) = memNoise(i,:)./std(memNoise(i,:));
        memNoise(i,:) = memNoise(i,:)*parameters.noiseScale(i);
    end
    
    % Make the necessary gamma pdfs for the current injected by spikes
    gammaMean = parameters.gammaMean./parameters.tau;
    gammaSTD = parameters.gammaSTD./parameters.tau;
    B = (gammaSTD(i)^2)/gammaMean(i);
    A = gammaMean(i)/B;
    y = gampdf(1:n,A,B); % mean A*B std sqrt(A*B^2)
    y = (y./max(y));
    maxL = find(y == max(y),1,'first');
    gammaFunL = maxL + find(y((maxL + 1):end) < 10^(-6),1,'first');
    gammaFun = y(1:gammaFunL);

    
    
    
    % Run the sensory neuron model
    [mV_E1,spkCt_E1,spkLoc_E1] = simpleModel_RS(I{1},T,parameters.tau,0,memNoise(1,:));
    
    
    % Figure out the scaling for the current from each spike based on the
    % decay
    B = (1/(1 - spkCt_E1))*log(parameters.minI/parameters.maxI);
    A = parameters.maxI*exp(B);
    weights = A*exp(-B*(1:spkCt_E1));
    
    % Add the gamma pdf pulses to the second neuron input current with the
    % decaying weights
    for iSpike = find((spkLoc_E1 + gammaFunL - 1) <= n)
        I{2}(spkLoc_E1(iSpike):(spkLoc_E1(iSpike) + gammaFunL - 1)) = ...
            I{2}(spkLoc_E1(iSpike):(spkLoc_E1(iSpike) + gammaFunL - 1)) + weights(iSpike)*gammaFun;
    end
    for iSpike = find((spkLoc_E1 + gammaFunL - 1) > n)
        I{2}(spkLoc_E1(iSpike):end) = ...
            I{2}(spkLoc_E1(iSpike):end) + weights(iSpike)*gammaFun(1:(n - spkLoc_E1(iSpike) + 1));
    end
    
    % Run the motor neuron model
    [mV_E2,spkCt_E2,spkLoc_E2] = simpleModel_RS(I{2},T,parameters.tau,0,memNoise(2,:));
    
    % Convert spike times to milliseconds
    spkLoc_E1 = spkLoc_E1*parameters.tau;
    spkLoc_E2 = spkLoc_E2*parameters.tau;
    
    % Record the results in the data structure
    data.spkTimes = cell([2,1]);
    data.spkTimes{1} = spkLoc_E1;
    data.spkTimes{2} = spkLoc_E2;
    data.parameters = parameters;
    data.weights = weights;
    data.currentRec = currentRec;
    data.currentOriginal = currentOriginal;
    data.stimTimes = stimTimes;
    data.stimTypes = stimTypes;
    data.time = time;
    data.expType = expType;
    
    
    
    
    
    
end





end