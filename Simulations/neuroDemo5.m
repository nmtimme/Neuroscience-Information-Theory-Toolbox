%% neuroDemo5
%
%   neuroDemo5 reproduces the canonical neuroscience experiment 
%   simulations. Note, the scratch directory in the simulations folder will
%   be  used to save data. These simulations will take up about 2.5 hours 
%   to run and it will require about 0.5 GB of space.
%
% Other m-files required: modeltype5, instinfo, inverseFNoise,
% simpleModels_FSI, simpleModel_RS
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% November 2016; Last revision: 1-Nov-2016

%% Set Parameters

disp('Note that this simulation will take about 2.5 hours to run and it requires about 0.5 GB of space.')

% Set the number of models to run for each subtype
nModelsPerSubtype = 20;

% Set the number of subtypes (this must be 6)
nSTs = 6;

% Set parameters for the subtypes
STparams = cell([nSTs,1]);
expType = cell([nSTs,1]);

% Center Out Task with Variable Response
parametersCO = struct;
parametersCO.maxFR = 100;
parametersCO.meanLead = -100;
parametersCO.stdLead = 10;
parametersCO.stdMP = 100;
parametersCO.nAngles = 8;
parametersCO.nNeurons = 100;
parametersCO.nTrials = 150;
parametersCO.intertrialT = 1000;
parametersCO.respRange = [0,1];
STparams{1} = parametersCO;
expType{1} = 'CenterOut';

% Center Out Task with Only Strong Responses
parametersCO = struct;
parametersCO.maxFR = 100;
parametersCO.meanLead = -100;
parametersCO.stdLead = 10;
parametersCO.stdMP = 100;
parametersCO.nAngles = 8;
parametersCO.nNeurons = 20;
parametersCO.nTrials = 150;
parametersCO.intertrialT = 1000;
parametersCO.respRange = [1,1];
STparams{2} = parametersCO;
expType{2} = 'CenterOut';

% Center Surround Cells
parametersCS = struct;
parametersCS.smallR = 0.1;
parametersCS.largeR = 0.3;
parametersCS.bgFR = 30;
parametersCS.centerFR = 100;
parametersCS.surroundFR = 1;
parametersCS.nStims = 400;
parametersCS.interstimT = 400;
parametersCS.intrastimT = 500;
parametersCS.nNeurons = 300;
STparams{3} = parametersCS;
expType{3} = 'CenterSurround';

% Place Cells with Only Strong Response
parametersPC = struct;
parametersPC.nNeurons = 200;
parametersPC.maxFR = 100;
parametersPC.bgFR = 20;
parametersPC.spatialSTD = 0.15;
parametersPC.stepDur = 10;
parametersPC.nSteps = 20000;
parametersPC.nStepsDim = 20;
parametersPC.respRange = [1,1];
STparams{4} = parametersPC;
expType{4} = 'PlaceCells';

% Sensory Habituation Model
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
STparams{5} = parametersSH;
expType{5} = 'SensoryHabit';

% Place Cells with Only Weak Response
parametersPC = struct;
parametersPC.nNeurons = 2;
parametersPC.maxFR = 100;
parametersPC.bgFR = 20;
parametersPC.spatialSTD = 0.15;
parametersPC.stepDur = 10;
parametersPC.nSteps = 20000;
parametersPC.nStepsDim = 20;
parametersPC.respRange = [0,0];
STparams{6} = parametersPC;
expType{6} = 'PlaceCells';

% Set the bin size in milliseconds for each subtype
BinSize = [25,25,25,100,50,100];

% Set the time window around the start of the stimuli to analyze in
% milliseconds
MaxLead = [600,600,200,0,200,0];
MaxLag = [400,400,700,BinSize(4),700,BinSize(6)];

% Set the number of 1D spatial bins for models that have spatial 
% distributions
nSBins = [NaN,NaN,4,4,NaN,4];

% Set the number of uniform count bins for each subtype
nUniCBins = [4,4,4,4,4,4];

% Get the scratch directory to use to store data
if ispc
    ScratchDir = [pwd,'\Scratch\ModelType5\'];
elseif isunix
    ScratchDir = [pwd,'/Scratch/ModelType5/'];
elseif ismac
    ScratchDir = [pwd,'/Scratch/ModelType5/'];
end

% Make the scratch directory if it does not exist
if exist(ScratchDir,'dir') ~= 7
    mkdir(ScratchDir)
end

% Find the file path for the analysis software and the simulation software
SimDir = pwd;
cd ..
AnaDir = pwd;
cd(SimDir)

%% Make the Model Data

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        data = modeltype5(expType{iSubType},'parameters', STparams{iSubType});
        
        save([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'.mat'],'data')
        
        disp(['Finished Generating Data for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Convert the Raw Data to DataRaster Format

% Go to the directory with the analysis software
cd(AnaDir)

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])
        
        % Find the number of neurons
        nNeurons = size(data.spkTimes,1);
        
        % Make the Data Raster
        DataRaster = cell([2,1]);
        
        if (iSubType == 1) || (iSubType == 2)
            
            % Center Out Movement Encoding
            
            % Format the Data
            [Temp,timeboundaries] = formattool(ones(size(data.spkTimes{1})),data.spkTimes{1},data.trialTimes,BinSize(iSubType),MaxLead(iSubType),MaxLag(iSubType),'count');
            DataRaster{1} = NaN([nNeurons,size(Temp,2),size(Temp,3)]);
            DataRaster{1}(1,:,:) = Temp;
            for i = 2:nNeurons
                DataRaster{1}(i,:,:) = formattool(ones(size(data.spkTimes{i})),data.spkTimes{i},data.trialTimes,BinSize(iSubType),MaxLead(iSubType),MaxLag(iSubType),'count');
            end
            DataRaster{2} = NaN([1,1,size(Temp,3)]);
            [waste,AngleStates] = histc(data.goalA,data.angles);
            DataRaster{2}(1,1,:) = reshape(AngleStates,[1,1,length(AngleStates)]);
            
            % Add parameters
            data.parameters.binsize = BinSize(iSubType);
            data.parameters.maxlead = MaxLead(iSubType);
            data.parameters.maxlag = MaxLag(iSubType);
            data.parameters.Time = timeboundaries;
            data.AngleStates = AngleStates;
            modelParams = struct;
            modelParams = data.parameters;
            
        elseif iSubType == 3
            
            % Center Surround Cell Data
            
            % Format the Data
            [Temp,timeboundaries] = formattool(ones(size(data.spkTimes{1})),data.spkTimes{1},data.stimTimes,BinSize(iSubType),MaxLead(iSubType),MaxLag(iSubType),'count');
            DataRaster{1} = NaN([nNeurons,size(Temp,2),size(Temp,3)]);
            DataRaster{1}(1,:,:) = Temp;
            for i = 2:nNeurons
                DataRaster{1}(i,:,:) = formattool(ones(size(data.spkTimes{i})),data.spkTimes{i},data.stimTimes,BinSize(iSubType),MaxLead(iSubType),MaxLag(iSubType),'count');
            end
            DataRaster{2} = NaN([1,1,size(Temp,3)]);
            spatialEdges = linspace(0,1,nSBins(iSubType) + 1);
            spatialEdges(end) = spatialEdges(end) + 0.1;
            [waste,xLocs] = histc(data.stimLocs(:,1),spatialEdges);
            [waste,yLocs] = histc(data.stimLocs(:,2),spatialEdges);
            [waste1,waste2,SpatialStates] = unique([xLocs,yLocs],'rows');
            DataRaster{2}(1,1,:) = reshape(SpatialStates,[1,1,length(SpatialStates)]);
            
            % Add parameters
            data.parameters.binsize = BinSize(iSubType);
            data.parameters.maxlead = MaxLead(iSubType);
            data.parameters.maxlag = MaxLag(iSubType);
            data.parameters.nSBins = nSBins(iSubType);
            data.parameters.Time = timeboundaries;
            data.parameters.cellLocs = data.cellLocs;
            data.parameters.stimLocs = data.stimLocs;
            data.parameters.SpatialStates = SpatialStates;
            modelParams = struct;
            modelParams = data.parameters;
            
        elseif (iSubType == 4) || (iSubType == 6)
            
            % Place Cell Data
            
            % Make a time bin vector
            binT = 1:BinSize(iSubType):size(data.animalLoc,1);
            
            % Format the Data
            [Temp,timeboundaries] = formattool(ones(size(data.spkTimes{1})),data.spkTimes{1},binT,BinSize(iSubType),MaxLead(iSubType),MaxLag(iSubType),'count');
            DataRaster{1} = NaN([nNeurons,size(Temp,2),size(Temp,3)]);
            DataRaster{1}(1,:,:) = Temp;
            for i = 2:nNeurons
                DataRaster{1}(i,:,:) = formattool(ones(size(data.spkTimes{i})),data.spkTimes{i},binT,BinSize(iSubType),MaxLead(iSubType),MaxLag(iSubType),'count');
            end
            DataRaster{2} = NaN([1,1,size(Temp,3)]);
            spatialEdges = linspace(0,1,nSBins(iSubType) + 1);
            spatialEdges(end) = spatialEdges(end) + 0.1;
            [waste,xLocs] = histc(data.animalLoc(:,1),spatialEdges);
            [waste,yLocs] = histc(data.animalLoc(:,2),spatialEdges);
            xLocs = squeeze(formattool(xLocs,1:length(xLocs),binT,BinSize(iSubType),MaxLead(iSubType),MaxLag(iSubType),'ave'));
            yLocs = squeeze(formattool(yLocs,1:length(yLocs),binT,BinSize(iSubType),MaxLead(iSubType),MaxLag(iSubType),'ave'));
            [waste1,waste2,SpatialStates] = unique(round([xLocs,yLocs]),'rows');
            DataRaster{2}(1,1,:) = reshape(SpatialStates,[1,1,length(SpatialStates)]);
            
            % Add parameters
            data.parameters.binsize = BinSize(iSubType);
            data.parameters.maxlead = MaxLead(iSubType);
            data.parameters.maxlag = MaxLag(iSubType);
            data.parameters.nSBins = nSBins(iSubType);
            data.parameters.Time = timeboundaries;
            data.parameters.binnedAnimalLocs = [xLocs,yLocs];
            data.parameters.cellPrefLocs = data.cellPrefLocs;
            modelParams = struct;
            modelParams = data.parameters;
            
            
        elseif iSubType == 5
            
            % Sensory Habituation Data
            
            % Format the Data
            [Temp,timeboundaries] = formattool(ones(size(data.spkTimes{1})),data.spkTimes{1},data.stimTimes,BinSize(iSubType),MaxLead(iSubType),MaxLag(iSubType),'count');
            DataRaster{1} = NaN([nNeurons,size(Temp,2),size(Temp,3)]);
            DataRaster{1}(1,:,:) = Temp;
            for i = 2:nNeurons
                DataRaster{1}(i,:,:) = formattool(ones(size(data.spkTimes{i})),data.spkTimes{i},data.stimTimes,BinSize(iSubType),MaxLead(iSubType),MaxLag(iSubType),'count');
            end
            DataRaster{2} = NaN([2,1,size(Temp,3)]);
            DataRaster{2}(1,1,:) = reshape(data.stimTypes,[1,1,length(data.stimTypes)]);
            DataRaster{2}(2,1,:) = reshape(1:length(data.stimTimes),[1,1,length(data.stimTimes)]);
            
            % Add parameters
            data.parameters.binsize = BinSize(iSubType);
            data.parameters.maxlead = MaxLead(iSubType);
            data.parameters.maxlag = MaxLag(iSubType);
            data.parameters.Time = timeboundaries;
            modelParams = struct;
            modelParams = data.parameters;
            
        end
        
        % Save the Data
        save([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'],'DataRaster','modelParams')
        
        disp(['Finished Converting Data to DataRaster Format for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Convert the DataRaster Format Data to StatesRasters

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])
        
        % Find the Number of neurons
        nNeurons = size(DataRaster{1},1);
        
        if (iSubType == 1) || (iSubType == 2)
            
            % Make the stating method assignment array
            MethodAssign = cell([nNeurons + 1,4]);
            MethodAssign{1,1} = 2;
            MethodAssign{1,2} = 1;
            MethodAssign{1,3} = 'Nat';
            MethodAssign{1,4} = {};
            for iNeuron = 1:nNeurons
                MethodAssign{iNeuron + 1,1} = 1;
                MethodAssign{iNeuron + 1,2} = iNeuron;
                MethodAssign{iNeuron + 1,3} = 'UniCB';
                MethodAssign{iNeuron + 1,4} = {nUniCBins(iSubType)};
            end
            
            modelParams.nUniCBins = nUniCBins(iSubType);
            
        elseif iSubType == 3
            
            % Make the stating method assignment array
            MethodAssign = cell([nNeurons + 1,4]);
            MethodAssign{1,1} = 2;
            MethodAssign{1,2} = 1;
            MethodAssign{1,3} = 'Nat';
            MethodAssign{1,4} = {};
            for iNeuron = 1:nNeurons
                MethodAssign{iNeuron + 1,1} = 1;
                MethodAssign{iNeuron + 1,2} = iNeuron;
                MethodAssign{iNeuron + 1,3} = 'UniCB';
                MethodAssign{iNeuron + 1,4} = {nUniCBins(iSubType)};
            end
            
            modelParams.nUniCBins = nUniCBins(iSubType);
            
        elseif (iSubType == 4) || (iSubType == 6)
            
            % Make the stating method assignment array
            MethodAssign = cell([nNeurons + 1,4]);
            MethodAssign{1,1} = 2;
            MethodAssign{1,2} = 1;
            MethodAssign{1,3} = 'Nat';
            MethodAssign{1,4} = {};
            for iNeuron = 1:nNeurons
                MethodAssign{iNeuron + 1,1} = 1;
                MethodAssign{iNeuron + 1,2} = iNeuron;
                MethodAssign{iNeuron + 1,3} = 'UniCB';
                MethodAssign{iNeuron + 1,4} = {nUniCBins(iSubType)};
            end
            
            modelParams.nUniCBins = nUniCBins(iSubType);
            
        elseif iSubType == 5
            
            % Make the stating method assignment array
            MethodAssign = cell([nNeurons + 2,4]);
            MethodAssign{1,1} = 2;
            MethodAssign{1,2} = 1;
            MethodAssign{1,3} = 'Nat';
            MethodAssign{1,4} = {};
            MethodAssign{2,1} = 2;
            MethodAssign{2,2} = 2;
            MethodAssign{2,3} = 'UniCB';
            MethodAssign{2,4} = {nUniCBins(iSubType)};
            for iNeuron = 1:nNeurons
                MethodAssign{iNeuron + 2,1} = 1;
                MethodAssign{iNeuron + 2,2} = iNeuron;
                MethodAssign{iNeuron + 2,3} = 'UniCB';
                MethodAssign{iNeuron + 2,4} = {nUniCBins(iSubType)};
            end
            
            modelParams.nUniCBins = nUniCBins(iSubType);
            
        end
        
        % Convert the data
        StatesRaster = data2states(DataRaster,MethodAssign);
        
        % Save the Data
        save([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'SR.mat'],'StatesRaster','modelParams')
        
        disp(['Finished Converting Data to StatesRaster Format for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Perform the Information Theory Analysis

InfoResults = cell([nSTs,nModelsPerSubtype]);
SubTypeParams = cell([nSTs,1]);
for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'SR.mat'])
        
        % Save the subtype parameters. We assume all runs of a subtype used the
        % same parameters.
        if iModel == 1
            SubTypeParams{iSubType} = modelParams;
        end
        
        % Figure out how many time steps we have in this data set
        nT = size(StatesRaster{1},2);
        
        % Figure out how many neurons we have in this data set
        nNeurons = size(StatesRaster{1},1);
        
        if iSubType == 1
            
            MIVals = NaN([nNeurons,nT]);
            
            Method = 'PairMI';
            for iT = 1:nT
                for iNeuron = 1:nNeurons
                    VariableIDs = {1,iNeuron,iT;2,1,1};
                    MIVals(iNeuron,iT) = instinfo(StatesRaster, Method, VariableIDs);
                end
            end
            
            InfoResults{iSubType,iModel} = {MIVals};
            
        elseif iSubType == 2
            
            PIDVals = NaN([nNeurons,nNeurons,nT,4]);
            
            Method = '2PID';
            for iT = 1:nT
                for iNeuron = 1:(nNeurons - 1)
                    for jNeuron = (iNeuron + 1):nNeurons
                        VariableIDs = {2,1,1;1,iNeuron,iT;1,jNeuron,iT};
                        PIDVals(iNeuron,jNeuron,iT,:) = instinfo(StatesRaster, Method, VariableIDs);
                    end
                end
            end
            
            InfoResults{iSubType,iModel} = {PIDVals};
            
        elseif iSubType == 3
            
            % Instead of pairing all neurons with the central neuron, draw random
            % pairs to avoid biases associated with spatial boundaries
            
            PairList = nchoosek(1:nNeurons,2);
            PairList = PairList(randi(size(PairList,1),[1,nNeurons]),:);
            
            PIDVals = NaN([nNeurons,nT,4]);
            
            Method = '2PID';
            for iT = 1:nT
                for iPair = 1:nNeurons
                    VariableIDs = {2,1,1;1,PairList(iPair,1),iT;1,PairList(iPair,2),iT};
                    PIDVals(iPair,iT,:) = instinfo(StatesRaster, Method, VariableIDs);
                end
            end
            
            MIVals = NaN([nNeurons,nT]);
            
            Method = 'PairMI';
            for iT = 1:nT
                for iNeuron = 1:nNeurons
                    VariableIDs = {1,iNeuron,iT;2,1,1};
                    MIVals(iNeuron,iT) = instinfo(StatesRaster, Method, VariableIDs);
                end
            end
            
            InfoResults{iSubType,iModel} = {PIDVals,MIVals,PairList};
            
        elseif iSubType == 4
            
            % Instead of pairing all neurons with the central neuron, draw random
            % pairs to avoid biases associated with spatial boundaries
            
            PairList = nchoosek(1:nNeurons,2);
            PairList = PairList(randi(size(PairList,1),[1,nNeurons]),:);
            
            PIDVals = NaN([nNeurons,4]);
            
            Method = '2PID';
            for iPair = 1:nNeurons
                VariableIDs = {2,1,1;1,PairList(iPair,1),1;1,PairList(iPair,2),1};
                PIDVals(iPair,:) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            MIVals = NaN([nNeurons,1]);
            
            Method = 'PairMI';
            for iNeuron = 1:nNeurons
                VariableIDs = {1,iNeuron,1;2,1,1};
                MIVals(iNeuron) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            InfoResults{iSubType,iModel} = {PIDVals,MIVals,PairList};
            
        elseif iSubType == 5
            
            MISE1Vals = NaN([1,nT]);
            MISE2Vals = NaN([1,nT]);
            MITE1Vals = NaN([1,nT]);
            MITE2Vals = NaN([1,nT]);
            
            Method = 'PairMI';
            for iT = 1:nT
                VariableIDs = {1,1,iT;2,1,1};
                MISE1Vals(iT) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            Method = 'PairMI';
            for iT = 1:nT
                VariableIDs = {1,2,iT;2,1,1};
                MISE2Vals(iT) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            Method = 'PairMI';
            for iT = 1:nT
                VariableIDs = {1,1,iT;2,2,1};
                MITE1Vals(iT) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            Method = 'PairMI';
            for iT = 1:nT
                VariableIDs = {1,2,iT;2,2,1};
                MITE2Vals(iT) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            InfoResults{iSubType,iModel} = {MISE1Vals,MISE2Vals,MITE1Vals,MITE2Vals};
            
        elseif iSubType == 6
            
            MIVals = NaN([nNeurons,1]);
            
            Method = 'PairMI';
            for iNeuron = 1:nNeurons
                VariableIDs = {1,iNeuron,1;2,1,1};
                MIVals(iNeuron) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            InfoResults{iSubType,iModel} = {MIVals};
            
        end
        
        disp(['Finished Analyzing for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

% Save the Results
save([ScratchDir,'ModelType5InfoResults.mat'],'InfoResults','SubTypeParams')

%% Gather Data for Plotting

% Get the Data for the center out model

% Find the file number for the first run of subtype 1, and load the data
iSubType = 1;
iModel = 1;
load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])

% First, get the data for the movement direction diagram
DP1 = struct;

xTemp = cos(data.angles);
yTemp = sin(data.angles);
DP1.points = [xTemp',yTemp'];
DP1.points = [0,0;DP1.points];

disp('Done with Data Packet 1')


% Second, get example firing rate profiles for a strong and weak encoders

% Find the example neurons
sNeuron = find(data.respS == max(data.respS));
wNeuron = find(data.respS == min(data.respS));

DP2 = struct;

% Load the data raster to get the firing rates
load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])

nT = size(DataRaster{1},2);
Directions = squeeze(DataRaster{2}(1,1,:));
Neurons = [sNeuron,wNeuron];

DP2.FRPlot = NaN([2,length(data.angles),nT,5]);
for iNeuron = 1:2
    for iDir = 1:length(data.angles)
        FRProfiles = squeeze(DataRaster{1}(Neurons(iNeuron),:,Directions == iDir))';
        FRProfiles = FRProfiles*(1000/modelParams.binsize);
        for iT = 1:nT
            temp = FRProfiles(:,iT);
            DP2.FRPlot(iNeuron,iDir,iT,1) = mean(temp);
            DP2.FRPlot(iNeuron,iDir,iT,2) = mean(temp) + std(temp);
            DP2.FRPlot(iNeuron,iDir,iT,3) = mean(temp) - std(temp);
            DP2.FRPlot(iNeuron,iDir,iT,4) = mean(temp) + (std(temp)/sqrt(length(temp)));
            DP2.FRPlot(iNeuron,iDir,iT,5) = mean(temp) - (std(temp)/sqrt(length(temp)));
        end
    end
end

% DP2.Time = mean(DataRaster.Time,2);
DP2.Time = mean(modelParams.Time,2);
DP2.neurons = Neurons;
DP2.respS = data.respS(Neurons);
DP2.prefA = data.prefA(Neurons);

disp('Done with Data Packet 2')



% Third, get the MI profiles for the example strong and weak encoders

% Load the information results
load([ScratchDir,'ModelType5InfoResults.mat'],'InfoResults')

DP3 = struct;

DP3.MIPlot = InfoResults{iSubType,iModel}{1}(Neurons,:);
DP3.Time = mean(modelParams.Time,2);
DP3.neurons = Neurons;
DP3.respS = data.respS(Neurons);
DP3.prefA = data.prefA(Neurons);

disp('Done with Data Packet 3')


% Fourth, get the MI values at max firing rate displacement as a function
% of response strength

% Load the first file to find parameters
load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])

% Figure out the total number of pairs of neurons
nNeurons = modelParams.nNeurons;
nNeuronsTotal = nModelsPerSubtype*nNeurons;

% Figure out the maximum point of the firing rate modification
temp = modelParams.Time(:,1) - modelParams.meanLead;
maxBin = find(temp <= 0,1,'last');

MIList = NaN([nNeuronsTotal,2]);

iNeuronTotal = 1;
for iModel = 1:nModelsPerSubtype
    load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])
    for iNeuron = 1:nNeurons
        MIList(iNeuronTotal,:) = [data.respS(iNeuron),InfoResults{iSubType,iModel}{1}(iNeuron,maxBin)];
        iNeuronTotal = iNeuronTotal + 1;
    end
end

% Make response strength edges
respEdges = 0:0.1:1;
[waste,Bin] = histc(MIList(:,1),respEdges);

DP4 = struct;

DP4.MICurves = NaN([length(respEdges) - 1,5]);
for iBin = 1:(length(respEdges) - 1)
    temp = MIList(Bin == iBin,2);
    DP4.MICurves(iBin,1) = mean(temp);
    DP4.MICurves(iBin,2) = mean(temp) + std(temp);
    DP4.MICurves(iBin,3) = mean(temp) - std(temp);
    DP4.MICurves(iBin,4) = mean(temp) + (std(temp)/sqrt(length(temp)));
    DP4.MICurves(iBin,5) = mean(temp) - (std(temp)/sqrt(length(temp)));
end
DP4.responses = respEdges(1:(end - 1)) + 0.5*(respEdges(2) - respEdges(1));

disp('Done with Data Packet 4')




% Fifth, get the PID values at max firing rate displacement as a function
% of angle in the strong response networks

% Find the subtype 2 files
iSubType = 2;

% Load the first file to find parameters
iModel = 1;
load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])

% Figure out the total number of pairs of neurons
nNeurons = modelParams.nNeurons;
nPairs = nModelsPerSubtype*nchoosek(nNeurons,2);

% Figure out the maximum point of the firing rate modification
temp = modelParams.Time(:,1) - modelParams.meanLead;
maxBin = find(temp <= 0,1,'last');

PIDList = NaN([nPairs,5]);

iPair = 1;
for iModel = 1:nModelsPerSubtype
    load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])
    for iNeuron = 1:(nNeurons - 1)
        for jNeuron = (iNeuron + 1):nNeurons
            temp = data.prefA(iNeuron) - (data.prefA(jNeuron) + [-2*pi,0,2*pi]);
            PIDList(iPair,:) = [temp(abs(temp) == min(abs(temp))),squeeze(InfoResults{iSubType,iModel}{1}(iNeuron,jNeuron,maxBin,:))'];
            iPair = iPair + 1;
        end
    end
end
PIDList(:,1) = PIDList(:,1)*(360/(2*pi));

% Make angle edges
angleEdges = -180:10:180;
[waste,Bin] = histc(PIDList(:,1),angleEdges);

DP5 = struct;

DP5.PIDCurves = NaN([length(angleEdges) - 1,4,5]);
for iBin = 1:(length(angleEdges) - 1)
    for iMeas = 1:4
        temp = PIDList(Bin == iBin,iMeas + 1);
        DP5.PIDCurves(iBin,iMeas,1) = mean(temp);
        DP5.PIDCurves(iBin,iMeas,2) = mean(temp) + std(temp);
        DP5.PIDCurves(iBin,iMeas,3) = mean(temp) - std(temp);
        DP5.PIDCurves(iBin,iMeas,4) = mean(temp) + (std(temp)/sqrt(length(temp)));
        DP5.PIDCurves(iBin,iMeas,5) = mean(temp) - (std(temp)/sqrt(length(temp)));
    end
end
DP5.angles = angleEdges(1:(end - 1)) + 0.5*(angleEdges(2) - angleEdges(1));


disp('Done with Data Packet 5')







% Get the data for the center surround model

% First, get data for the diagram and the example stimulations

% Load the data
iSubType = 3;
iModel = 1;
load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])

% Calculate the radius of all the stimuli from the first neuron
r = sqrt(sum((data.stimLocs - repmat([0.5,0.5],[data.parameters.nStims,1])).^2,2));

% Find a stimulus in the center, in the middle of the surround, and in the
% middle of the exterior
specStim = zeros([3,1]);
specStim(1) = find(r == min(r),1,'first');
specStim(2) = find(abs(r - mean([data.parameters.smallR,data.parameters.largeR])) == min(abs(r - mean([data.parameters.smallR,data.parameters.largeR]))),1,'first');
specStim(3) = find(abs(r - mean([0.5,data.parameters.largeR])) == min(abs(r - mean([0.5,data.parameters.largeR]))),1,'first');

% Get spikes for these three special cases
spikes = cell([3,1]);
for iCase = 1:3
    temp = data.spkTimes{1};
    temp(temp > (data.stimTimes(specStim(iCase)) + data.parameters.intrastimT + 0.5*data.parameters.interstimT)) = [];
    temp(temp < (data.stimTimes(specStim(iCase)) - 0.5*data.parameters.interstimT)) = [];
    spikes{iCase} = temp - data.stimTimes(specStim(iCase));
end

DP6 = struct;

DP6.radii = [data.parameters.smallR,data.parameters.largeR];
DP6.stimLocs = data.stimLocs(specStim,:);
DP6.timeLocs = [-0.5*data.parameters.interstimT,0,data.parameters.intrastimT,data.parameters.intrastimT + 0.5*data.parameters.interstimT];
DP6.spikes = spikes;

disp('Done with Data Packet 6')


% Get the information for the box plot of the stimulus location encoding
iSubType = 3;
iModel = 1;
load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])

% Find the number of time bins
nT = size(InfoResults{iSubType,iModel}{2},2);

MIValues = NaN([nModelsPerSubtype,nT]);
for iModel = 1:nModelsPerSubtype
    MIValues(iModel,:) = InfoResults{iSubType,iModel}{2}(1,:);
end

errorBarVals = NaN([nT,5]);
for iT = 1:nT
    temp = MIValues(:,iT);
    errorBarVals(iT,1) = mean(temp);
    errorBarVals(iT,2) = mean(temp) + std(temp);
    errorBarVals(iT,3) = mean(temp) - std(temp);
    errorBarVals(iT,4) = mean(temp) + (std(temp)/sqrt(length(temp)));
    errorBarVals(iT,5) = mean(temp) - (std(temp)/sqrt(length(temp)));
end

DP7 = struct;
DP7.errorBarVals = errorBarVals;
DP7.xLims = [-modelParams.maxlead,modelParams.maxlag];
DP7.time = mean(modelParams.Time,2);

disp('Done with Data Packet 7')




% Get the PID value heat maps

% Find the time bins during the stimulus
goodTBins = find((mean(modelParams.Time,2) > 0) & (mean(modelParams.Time,2) < modelParams.intrastimT));

% Set the number of spatial bins for the heat maps and make the edges
% vector
nSBins = 20;
sEdges = linspace(0,1,nSBins + 1);
sEdges(end) = sEdges(end) + 0.1;

PIDHeatMaps = zeros([nSBins,nSBins,4]);
PIDCounts = zeros([nSBins,nSBins,4]);
for iModel = 1:nModelsPerSubtype
    
    load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])
    
    % Go through each neuron and add the data to the heat maps
    for iNeuron = 1:modelParams.nNeurons
        
        % Get both neurons' locations
        xFocusOrig = modelParams.cellLocs(InfoResults{iSubType,iModel}{3}(iNeuron,1),1);
        yFocusOrig = modelParams.cellLocs(InfoResults{iSubType,iModel}{3}(iNeuron,1),2);
        xVarOrig = modelParams.cellLocs(InfoResults{iSubType,iModel}{3}(iNeuron,2),1);
        yVarOrig = modelParams.cellLocs(InfoResults{iSubType,iModel}{3}(iNeuron,2),2);
        
        % Transform the coordinates so the focus neuron is at (0.5,0.5)
        xVarShift = xVarOrig + (0.5 - xFocusOrig);
        yVarShift = yVarOrig + (0.5 - yFocusOrig);
        
        % Transform the coordinates back into the unit square, if necessary
        if xVarShift > 1
            xVarShift = xVarShift - 1;
        elseif xVarShift < 0
            xVarShift = xVarShift + 1;
        end
        if yVarShift > 1
            yVarShift = yVarShift - 1;
        elseif yVarShift < 0
            yVarShift = yVarShift + 1;
        end
        
        % Calculate the spatial bin for the variable neuron
        [waste,xBin] = histc(xVarShift,sEdges);
        [waste,yBin] = histc(yVarShift,sEdges);
        
        for iT = 1:length(goodTBins)
            for iMeas = 1:4
                PIDHeatMaps(xBin,yBin,iMeas) = PIDHeatMaps(xBin,yBin,iMeas) + InfoResults{iSubType,iModel}{1}(iNeuron,goodTBins(iT),iMeas);
                PIDCounts(xBin,yBin,iMeas) = PIDCounts(xBin,yBin,iMeas) + 1;
            end
        end
    end
end

DP8 = struct;
DP8.PIDHeatMaps = PIDHeatMaps;
DP8.PIDCounts = PIDCounts;

disp('Done with Data Packet 8')



% Make PID plots as a function of radial distance

% Set the size of the spatial bins
sBins = 0.04;
sEdges = 0:sBins:2;

PIDRadVals = NaN([300*length(goodTBins)*20,length(sEdges),4]);
PIDRadCounts = ones([length(sEdges),1]);
for iModel = 1:nModelsPerSubtype
    
    load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])
    
    for iNeuron = 1:modelParams.nNeurons
    
        % Get both neurons' locations
        xFocusOrig = modelParams.cellLocs(InfoResults{iSubType,iModel}{3}(iNeuron,1),1);
        yFocusOrig = modelParams.cellLocs(InfoResults{iSubType,iModel}{3}(iNeuron,1),2);
        xVarOrig = modelParams.cellLocs(InfoResults{iSubType,iModel}{3}(iNeuron,2),1);
        yVarOrig = modelParams.cellLocs(InfoResults{iSubType,iModel}{3}(iNeuron,2),2);
        
        % Figure out the minimum distance given the periodic boundary
        % conditions
        xVarTile = xVarOrig + [-1,-1,-1,0,0,0,1,1,1];
        yVarTile = yVarOrig + [-1,0,1,-1,0,1,-1,0,1];
        r = min(sqrt((xFocusOrig*ones([1,9]) - xVarTile).^2 + (yFocusOrig*ones([1,9]) - yVarTile).^2));
        
        % Calculate the spatial bin for the variable neuron
        [waste,rBin] = histc(r,sEdges);
        
        for iT = 1:length(goodTBins)
            for iMeas = 1:4
                PIDRadVals(PIDRadCounts(rBin),rBin,iMeas) = InfoResults{iSubType,iModel}{1}(iNeuron,goodTBins(iT),iMeas);
            end
            PIDRadCounts(rBin) = PIDRadCounts(rBin) + 1;
        end
        
    end
    
end

PIDRadVals(:,PIDRadCounts == 1,:) = [];
sEdges(PIDRadCounts == 1) = [];

DP15 = struct;

DP15.PIDCurves = NaN([length(sEdges),4,5]);
for iBin = 1:length(sEdges)
    for iMeas = 1:4
        temp = PIDRadVals(:,iBin,iMeas);
        temp(~isfinite(temp)) = [];
        DP15.PIDCurves(iBin,iMeas,1) = mean(temp);
        DP15.PIDCurves(iBin,iMeas,2) = mean(temp) + std(temp);
        DP15.PIDCurves(iBin,iMeas,3) = mean(temp) - std(temp);
        DP15.PIDCurves(iBin,iMeas,4) = mean(temp) + (std(temp)/sqrt(length(temp)));
        DP15.PIDCurves(iBin,iMeas,5) = mean(temp) - (std(temp)/sqrt(length(temp)));
    end
end
DP15.dist = sEdges + 0.5*sBins;

disp('Done with Data Packet 15')








% Gather the data for the place cell figures

% First, make a spike rate heat map and a linger time heat map for an
% example neuron (the first neuron because we forced it to be in the
% center)

% Load the data
iSubType = 4;
iModel = 1;
load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])

% Set the number of spatial bins
nSBins = 20;

% Make the spatial edges vector
sptEdges = linspace(0,1,nSBins + 1);
sptEdges(end) = sptEdges(end) + 0.1;

% Find the locations of the animal during spikes
spkLocs = data.animalLoc(data.spkTimes{1},:);

% Make histograms
[waste,xSpks] = histc(spkLocs(:,1),sptEdges);
[waste,ySpks] = histc(spkLocs(:,2),sptEdges);
[waste,xAn] = histc(data.animalLoc(:,1),sptEdges);
[waste,yAn] = histc(data.animalLoc(:,2),sptEdges);
spkHM = zeros([nSBins,nSBins]);
anHM = zeros([nSBins,nSBins]);
for i = 1:length(xSpks)
    spkHM(xSpks(i),ySpks(i)) = spkHM(xSpks(i),ySpks(i)) + 1;
end
for i = 1:length(xAn)
    anHM(xAn(i),yAn(i)) = anHM(xAn(i),yAn(i)) + 1;
end

% Convert linger time to seconds from milliseconds
anHM = anHM*0.001;

% Convert to spike rate
spkHM = spkHM./anHM;

% Put the results in a structure
DP9 = struct;
DP9.spkHM = spkHM;
DP9.anHM = anHM;

disp('Done with Data Packet 9')




% Make a box plot of the position encoding by the center neuron across data
% sets

MIValues = NaN([nModelsPerSubtype,2]);

iSubType = 4;
for iModel = 1:nModelsPerSubtype
    MIValues(iModel,1) = InfoResults{iSubType,iModel}{2}(1);
end
iSubType = 6;
for iModel = 1:nModelsPerSubtype
    MIValues(iModel,2) = InfoResults{iSubType,iModel}{1}(1);
end

errorBarVals = NaN([5,2]);
tempMIVals = cell([1,2]);
for iType = 1:2
    temp = sort(MIValues(:,iType),'ascend');
    errorBarVals(1,iType) = temp(1);
    errorBarVals(2,iType) = temp(ceil(0.25*length(temp)));
    errorBarVals(3,iType) = median(temp);
    errorBarVals(4,iType) = temp(floor(0.75*length(temp)));
    errorBarVals(5,iType) = temp(end);
    tempMIVals{iType} = temp;
end
pTemp = ranksum(tempMIVals{1},tempMIVals{2});

DP10 = struct;
DP10.errorBarVals = errorBarVals;
DP10.sigDots = 0;
if pTemp < 0.05
    DP10.sigDots = DP10.sigDots + 1;
end
if pTemp < 0.01
    DP10.sigDots = DP10.sigDots + 1;
end
if pTemp < 0.001
    DP10.sigDots = DP10.sigDots + 1;
end

disp('Done with Data Packet 10')




% Make heat maps of the PID values between the center neuron and the others

iSubType = 4;

% Set the number of spatial bins for the heat maps and make the edges
% vector
nSBins = 20;
sEdges = linspace(0,1,nSBins + 1);
sEdges(end) = sEdges(end) + 0.1;

PIDHeatMaps = zeros([nSBins,nSBins,4]);
PIDCounts = zeros([nSBins,nSBins,4]);
for iModel = 1:nModelsPerSubtype
    load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'],'modelParams')
    
    for iNeuron = 1:modelParams.nNeurons
    
        % Get both neurons' locations
        xFocusOrig = modelParams.cellPrefLocs(InfoResults{iSubType,iModel}{3}(iNeuron,1),1);
        yFocusOrig = modelParams.cellPrefLocs(InfoResults{iSubType,iModel}{3}(iNeuron,1),2);
        xVarOrig = modelParams.cellPrefLocs(InfoResults{iSubType,iModel}{3}(iNeuron,2),1);
        yVarOrig = modelParams.cellPrefLocs(InfoResults{iSubType,iModel}{3}(iNeuron,2),2);
        
        % Transform the coordinates so the focus neuron is at (0.5,0.5)
        xVarShift = xVarOrig + (0.5 - xFocusOrig);
        yVarShift = yVarOrig + (0.5 - yFocusOrig);
        
        % Transform the coordinates back into the unit square, if necessary
        if xVarShift > 1
            xVarShift = xVarShift - 1;
        elseif xVarShift < 0
            xVarShift = xVarShift + 1;
        end
        if yVarShift > 1
            yVarShift = yVarShift - 1;
        elseif yVarShift < 0
            yVarShift = yVarShift + 1;
        end
        
        % Calculate the spatial bin for the variable neuron
        [waste,xBin] = histc(xVarShift,sEdges);
        [waste,yBin] = histc(yVarShift,sEdges);
        
        for iMeas = 1:4
            PIDHeatMaps(xBin,yBin,iMeas) = PIDHeatMaps(xBin,yBin,iMeas) + InfoResults{iSubType,iModel}{1}(iNeuron,iMeas);
            PIDCounts(xBin,yBin,iMeas) = PIDCounts(xBin,yBin,iMeas) + 1;
        end
    end
    
end

DP11 = struct;
DP11.PIDHeatMaps = PIDHeatMaps./PIDCounts;

disp('Done with Data Packet 11')


% Make PID plots as a function of radial distance



% Set the size of the spatial bins
sBins = 0.04;
sEdges = 0:sBins:2;

PIDRadVals = NaN([200*20,length(sEdges),4]);
PIDRadCounts = ones([length(sEdges),1]);
for iModel = 1:nModelsPerSubtype
    load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'],'modelParams')
    
    for iNeuron = 1:modelParams.nNeurons
    
        % Get both neurons' locations
        xFocusOrig = modelParams.cellPrefLocs(InfoResults{iSubType,iModel}{3}(iNeuron,1),1);
        yFocusOrig = modelParams.cellPrefLocs(InfoResults{iSubType,iModel}{3}(iNeuron,1),2);
        xVarOrig = modelParams.cellPrefLocs(InfoResults{iSubType,iModel}{3}(iNeuron,2),1);
        yVarOrig = modelParams.cellPrefLocs(InfoResults{iSubType,iModel}{3}(iNeuron,2),2);
        
        % Figure out the minimum distance given the periodic boundary
        % conditions
        xVarTile = xVarOrig + [-1,-1,-1,0,0,0,1,1,1];
        yVarTile = yVarOrig + [-1,0,1,-1,0,1,-1,0,1];
        r = min(sqrt((xFocusOrig*ones([1,9]) - xVarTile).^2 + (yFocusOrig*ones([1,9]) - yVarTile).^2));
        
        % Calculate the spatial bin for the variable neuron
        [waste,rBin] = histc(r,sEdges);
        
        for iMeas = 1:4
            PIDRadVals(PIDRadCounts(rBin),rBin,iMeas) = InfoResults{iSubType,iModel}{1}(iNeuron,iMeas);
        end
        PIDRadCounts(rBin) = PIDRadCounts(rBin) + 1;
    end
    
end

PIDRadVals(:,PIDRadCounts == 1,:) = [];
sEdges(PIDRadCounts == 1) = [];

DP14 = struct;

DP14.PIDCurves = NaN([length(sEdges),4,5]);
for iBin = 1:length(sEdges)
    for iMeas = 1:4
        temp = PIDRadVals(:,iBin,iMeas);
        temp(~isfinite(temp)) = [];
        DP14.PIDCurves(iBin,iMeas,1) = mean(temp);
        DP14.PIDCurves(iBin,iMeas,2) = mean(temp) + std(temp);
        DP14.PIDCurves(iBin,iMeas,3) = mean(temp) - std(temp);
        DP14.PIDCurves(iBin,iMeas,4) = mean(temp) + (std(temp)/sqrt(length(temp)));
        DP14.PIDCurves(iBin,iMeas,5) = mean(temp) - (std(temp)/sqrt(length(temp)));
    end
end
DP14.dist = sEdges + 0.5*sBins;

disp('Done with Data Packet 14')







% Gather results for the sensory habituation model

% Find the file number for the first run of subtype 5, and load the data
iSubType = 5;
iModel = 1;
load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])

% Figure out which stimulations to look at 
specStims = zeros([1,3]);
specStims(1) = 1;
specStims(2) = round(length(data.stimTypes)/2);
while data.stimTypes(specStims(2)) ~= 1
    specStims(2) = specStims(2) + 1;
end
specStims(3) = find(data.stimTypes == 1,1,'last');

% Get spikes for these three special cases
spikes = cell([3,2]);
for iCase = 1:3
    for iNeuron = 1:2
        temp = data.spkTimes{iNeuron};
        temp(temp > (data.stimTimes(specStims(iCase)) + 2*data.parameters.intrastimT + 1.5*data.parameters.interstimT)) = [];
        temp(temp < (data.stimTimes(specStims(iCase)) - 0.5*data.parameters.interstimT)) = [];
        spikes{iCase,iNeuron} = temp - data.stimTimes(specStims(iCase));
    end
end

DP12 = struct;

DP12.timeLocs = [-0.5*data.parameters.interstimT,...
    0,...
    data.parameters.intrastimT,...
    data.parameters.intrastimT + data.parameters.interstimT,...
    2*data.parameters.intrastimT + data.parameters.interstimT,...
    2*data.parameters.intrastimT + 1.5*data.parameters.interstimT];
DP12.spikes = spikes;

disp('Done with Data Packet 12')



% Get the information for the box plot of the stimulus and trial number
% encoding

% Find the number of time bins
nT = size(InfoResults{iSubType,iModel}{2},2);

MIValues = NaN([nModelsPerSubtype,nT,4]);
for iModel = 1:nModelsPerSubtype
    MIValues(iModel,:,1) = InfoResults{iSubType,iModel}{1}; % Stimulus E1 encoding
    MIValues(iModel,:,2) = InfoResults{iSubType,iModel}{2}; % Stimulus E2 encoding
    MIValues(iModel,:,3) = InfoResults{iSubType,iModel}{3}; % Trail number E1 encoding
    MIValues(iModel,:,4) = InfoResults{iSubType,iModel}{4}; % Trial number E2 encoding
end

errorBarVals = NaN([nT,4,5]);
for iT = 1:nT
    for iMeas = 1:4
        temp = MIValues(:,iT,iMeas);
        errorBarVals(iT,iMeas,1) = mean(temp);
        errorBarVals(iT,iMeas,2) = mean(temp) + std(temp);
        errorBarVals(iT,iMeas,3) = mean(temp) - std(temp);
        errorBarVals(iT,iMeas,4) = mean(temp) + (std(temp)/sqrt(length(temp)));
        errorBarVals(iT,iMeas,5) = mean(temp) - (std(temp)/sqrt(length(temp)));
    end
end

DP13 = struct;
DP13.errorBarVals = errorBarVals;
DP13.xLims = [-modelParams.maxlead,modelParams.maxlag];
load([ScratchDir,'ModelType5ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])
DP13.time = mean(modelParams.Time,2);

disp('Done with Data Packet 13')


% Save the Results
save([ScratchDir,'ModelType5PlotData'],'DP1','DP2','DP3','DP4','DP5','DP6','DP7','DP8','DP9','DP10','DP11','DP12','DP13','DP14','DP15')


%% Make the Center-Out Task Figure

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.18;

% Set the top and bottom margins in inches
tmargin = 0.4;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace1 = 0.4;
hspace2 = 0.15;

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set the height and width of the direction diagram (square)
width1 = 1.75;

% Set Paper size
papersize = [7.4 6];

% calculate the width and height based on dimensions above in inches
width2 = (papersize(1) - lmargin - rmargin - width1 - hspace1 - hspace2)/2;
width3 = (papersize(1) - lmargin - rmargin - 2*hspace1)/3;
width4 = (papersize(1) - lmargin - rmargin - 3*hspace1)/4;
height1 = width1;
height2 = (papersize(2) - tmargin - bmargin - 2*vspace - height1)/2;

LeftCoord1 = [lmargin,lmargin + width1 + hspace1,lmargin + width1 + width2 + hspace1 + hspace2];
LeftCoord2 = (0:(3 - 1))*(width3 + hspace1) + lmargin;
LeftCoord3 = (0:(4 - 1))*(width4 + hspace1) + lmargin;
BottomCoord = [papersize(2) - tmargin - height1,bmargin + vspace + height2,bmargin];

% Convert to fraction of page sizes
LeftCoord1 = LeftCoord1/papersize(1);
LeftCoord2 = LeftCoord2/papersize(1);
LeftCoord3 = LeftCoord3/papersize(1);
BottomCoord = BottomCoord/papersize(2);
width1 = width1/papersize(1);
width2 = width2/papersize(1);
width3 = width3/papersize(1);
width4 = width4/papersize(1);
height1 = height1/papersize(2);
height2 = height2/papersize(2);

% Set Colors
PColor = {[0,0,1],[1,0,0],[0,210/255,50/255],[150/255,0,200/255],[238/255,130/255,238/255],[200/255,0,255/255],[0/255,255/255,185/255],[238/255,130/255,238/255]};
PColorGS = cell([5,1]);
for i = 1:5
    PColorGS{i} = [0,0,0] + 0.15*(i - 1);
end

% Make the colors for the directions
dirColors = jet(size(DP1.points,1) - 1);

% Set the background gray color
BGColor = 0.85*ones([1,3]);

% Set fontsizes
LegFS = 6;
TitleFS = 8;
AxLabelFS = 7;
UnitFS = 6;
StimLabelFS = 6;
AngleNameFS = 8;

% Set the linewidths
VoltageLW = 0.6;
BPThickLW = 3;
BPThinLW = 0.6;
StimLineLW = 0.8;
DelayLW = 2;
AxLW = 0.4;
SpikeLW = 0.4;
FRLW = 1;
DirLW = 5;
MIRespLW = 0.8;

% Set dot sizes
BPDotSz = 8;
DirDotSz = 200;
MIRespDotSz = 10;

% Make PID names
PIDNamesLong = {'Redundancy','Unique A','Unique B','Synergy'};
PIDNamesShort = {'Red (bits)','Uniq A (bits)','Uniq B (bits)','Syn (bits)'};


% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');



% Make the direction diagram

% Set the zoom factors
aZoom = 0.4;

f = subplot('Position',[LeftCoord1(1),BottomCoord(1),width1,height1]);
hold on

DP1.points = (aZoom*DP1.points) + 0.5;

for iDir = 1:(size(DP1.points,1) - 1)
    line([DP1.points(1,1),DP1.points(iDir + 1,1)],[DP1.points(1,2),DP1.points(iDir + 1,2)],'Color',dirColors(iDir,:),'LineWidth',DirLW)
end
scatter(DP1.points(1,1),DP1.points(1,2),DirDotSz,'k','filled')

xlim([0,1])
ylim([0,1])

title('Directions of Motion','FontSize',TitleFS)
axis off



% Plot the firing rate profiles

yMax = max(max(max(squeeze(DP2.FRPlot(:,:,:,1)))));
yMin = min(min(min(squeeze(DP2.FRPlot(:,:,:,1)))));
dy = yMax - yMin;

for iNeuron = 1:2
    
    f = subplot('Position',[LeftCoord1(1 + iNeuron),BottomCoord(1),width2,height1]);
    hold on
    
    line([0,0],[yMin - 0.1*dy,yMax + 0.1*dy],'Color','k','LineStyle','--','LineWidth',StimLineLW)
    text(0 + 0.02*(DP2.Time(end) - DP2.Time(1)),yMax,'Movement Onset','FontSize',LegFS)
    
    for iDir = 1:size(DP2.FRPlot,2)
        plot(DP2.Time,squeeze(DP2.FRPlot(iNeuron,iDir,:,1)),'Color',dirColors(iDir,:),'LineWidth',FRLW)
    end
    
    xlim([DP2.Time(1),DP2.Time(end)])
    ylim([yMin - 0.1*dy,yMax + 0.1*dy])
    
    set(gca,'FontSize',UnitFS)
    set(gca,'Layer','Top')
    
    xlabel('Time (ms)','FontSize',AxLabelFS)
    
    if iNeuron == 1
        ylabel('Firing Rate (Hz)','FontSize',AxLabelFS)
        title(['Example Strong Encoder Firing Rate',char(10),'r = ',num2str(DP2.respS(iNeuron),2),', \theta_{pref} = ',num2str(DP2.prefA(iNeuron)*360/(2*pi),3),char(176)],'FontSize',TitleFS)
    else
        title(['Example Weak Encoder Firing Rate',char(10),'r = ',num2str(DP2.respS(iNeuron),2),', \theta_{pref} = ',num2str(DP2.prefA(iNeuron)*360/(2*pi),3),char(176)],'FontSize',TitleFS)
        set(gca,'YTickLabel',{})
    end
    
end


% Plot the max FR deviation MI values across response

f = subplot('Position',[LeftCoord1(1),BottomCoord(2),width1,height2]);
hold on

yMax = max(max(DP4.MICurves(:,1:3)));
yMin = min(min(DP4.MICurves(:,1:3)));
dy = yMax - yMin;

scatter(DP4.responses,DP4.MICurves(:,1),MIRespDotSz,'k','filled')
for iResp = 1:length(DP4.responses)
    line([DP4.responses(iResp),DP4.responses(iResp)],[DP4.MICurves(iResp,2),DP4.MICurves(iResp,3)],'Color','k','LineWidth',MIRespLW)
end

xlim([0,1])
ylim([yMin - 0.1*dy,yMax + 0.1*dy])

set(gca,'FontSize',UnitFS)

xlabel('Response r','FontSize',AxLabelFS)
ylabel('Mutual Information (bits)','FontSize',AxLabelFS)
title('Encoding Response Dependency','FontSize',TitleFS)



% Plot the MI curves for the examples

yMax = max(DP3.MIPlot(:));
yMin = min(DP3.MIPlot(:));
dy = yMax - yMin;

for iNeuron = 1:2
    
    f = subplot('Position',[LeftCoord1(1 + iNeuron),BottomCoord(2),width2,height2]);
    hold on
    
    line([0,0],[yMin - 0.1*dy,yMax + 0.1*dy],'Color','k','LineStyle','--','LineWidth',StimLineLW)
    text(0 + 0.02*(DP2.Time(end) - DP2.Time(1)),yMax,'Movement Onset','FontSize',LegFS)
    
    plot(DP3.Time,squeeze(DP3.MIPlot(iNeuron,:)),'Color','k','LineWidth',FRLW)
    
    xlim([DP2.Time(1),DP2.Time(end)])
    ylim([yMin - 0.1*dy,yMax + 0.1*dy])
    
    set(gca,'FontSize',UnitFS)
    set(gca,'Layer','Top')
    
    xlabel('Time (ms)','FontSize',AxLabelFS)
    
    if iNeuron == 1
        ylabel('Mutual Information (Bits)','FontSize',AxLabelFS)
        title('Example Strong Encoder MI','FontSize',TitleFS)
    else
        title('Example Weak Encoder MI','FontSize',TitleFS)
        set(gca,'YTickLabel',{})
    end
    
end

% Plot the PID curves

yMax = max(max(max(DP5.PIDCurves(:,:,1:3))));
yMin = min(min(min(DP5.PIDCurves(:,:,1:3))));
dy = yMax - yMin;

for iMeas = 1:4
    
    f = subplot('Position',[LeftCoord3(iMeas),BottomCoord(3),width4,height2]);
    hold on
    
    line([0,0],[yMin - 0.1*dy,yMax + 0.1*dy],'Color','k','LineStyle','--','LineWidth',StimLineLW)
    
    DP5.angles = DP5.angles;
    
    scatter(DP5.angles,DP5.PIDCurves(:,iMeas,1),MIRespDotSz,'k','filled')
    for iAng = 1:length(DP5.angles)
        line([DP5.angles(iAng),DP5.angles(iAng)],[DP5.PIDCurves(iAng,iMeas,2),DP5.PIDCurves(iAng,iMeas,3)],'Color','k','LineWidth',MIRespLW)
    end
    
    xlim([-180,180])
    ylim([yMin - 0.1*dy,yMax + 0.1*dy])
    
    set(gca,'FontSize',UnitFS)
    set(gca,'Layer','Top')
    
    xlabel('Relative Angle (degrees)','FontSize',AxLabelFS)
    ylabel(PIDNamesShort{iMeas},'FontSize',AxLabelFS)
    title(PIDNamesLong{iMeas},'FontSize',TitleFS)
    
end


% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo5CO')
cd(SimDir)





%% Make the Center Surround Retinal Ganglion Cell Figure

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.18;

% Set the top and bottom margins in inches
tmargin = 0.4;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace1 = 0.2;
hspace2 = 0.45;
hspace3 = 0.5;

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set the height and width of the direction diagram (square) in inches
width1 = 1.75;

% Set the width of the colorbar in inches
widthCB = 0.2;

% Set Paper size
papersize = [7.4 6];

% calculate the width and height based on dimensions above in inches
width2 = (papersize(1) - lmargin - rmargin - width1 - hspace1 - hspace2)/2;
width3 = (papersize(1) - lmargin - 3*hspace3 - rmargin)/4;
height1 = width1;
height2 = width3;

LeftCoord1 = [lmargin,lmargin + width1 + hspace1,lmargin + width1 + width2 + hspace1 + hspace2];
LeftCoord2 = (0:(4 - 1))*(width3 + hspace3) + lmargin;
BottomCoord = [papersize(2) - tmargin - height1,papersize(2) - tmargin - height1 - vspace - height2];

% Convert to fraction of page sizes
LeftCoord1 = LeftCoord1/papersize(1);
LeftCoord2 = LeftCoord2/papersize(1);
BottomCoord = BottomCoord/papersize(2);
width1 = width1/papersize(1);
width2 = width2/papersize(1);
width3 = width3/papersize(1);
widthCB = widthCB/papersize(1);
height1 = height1/papersize(2);
height2 = height2/papersize(2);

% Set Colors
PColor = {[0,0,1],[1,0,0],[0,210/255,50/255],[150/255,0,200/255],[238/255,130/255,238/255],[200/255,0,255/255],[0/255,255/255,185/255],[238/255,130/255,238/255]};
PColorGS = cell([5,1]);
for i = 1:5
    PColorGS{i} = [0,0,0] + 0.15*(i - 1);
end

% Set the background gray color
BGColor = 0.85*ones([1,3]);

% Set fontsizes
LegFS = 6;
TitleFS = 8;
AxLabelFS = 7;
UnitFS = 6;
StimLabelFS = 8;
AngleNameFS = 8;
StimDotFS = 8;
RecLabelFS = 8;
CBFS = 6;

% Set the linewidths
VoltageLW = 0.6;
BPThickLW = 3;
BPThinLW = 0.6;
StimLineLW = 0.8;
DelayLW = 2;
AxLW = 0.4;
SpikeLW = 0.4;
FRLW = 1;
DirLW = 5;
MIRespLW = 0.8;
CircleLW = 1.2;
ErrBarLW = 0.8;
PIDLW = 1;

% Set dot sizes
StimDotSz = 20;
MIDotSz = 15;
PIDDotSz = 20;

% Make PID names
PIDNamesLong = {'Redundancy','Unique Center','Unique Variable','Synergy'};
PIDNamesShort = {'Red (bits)','Uniq A (bits)','Uniq B (bits)','Syn (bits)'};

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');



% Make the center surround diagram

f = subplot('Position',[LeftCoord1(1),BottomCoord(1),width1,height1]);
hold on

% Figure out the point size of the two dots
DotSz = zeros([1,2]);
for i = 1:2
    a = (2*DP6.radii(i))^2; % area of the box in units of the 1 by 1 space
    a = a*width1*papersize(1)*height1*papersize(2); % area in units of inches (because the space is 1 by 1)
    a = a/(0.0138889^2); % area in points (1 point = 0.0138889 inches)
    DotSz(i) = a;
end

% Plot the region boundaries
scatter(0.5,0.5,DotSz(1),[0.6,0.6,0.6],'LineWidth',CircleLW)
scatter(0.5,0.5,DotSz(2),[0.6,0.6,0.6],'LineWidth',CircleLW)

% text(0.5,0.5,'+','FontSize',RecLabelFS,'FontSize',RecLabelFS,'HorizontalAlignment','Center','VerticalAlignment','Middle')
text(0.46,0.46,'+','FontSize',RecLabelFS,'FontSize',RecLabelFS,'HorizontalAlignment','Center','VerticalAlignment','Middle')
text(0.5 - mean(DP6.radii)*cosd(45),0.5 + mean(DP6.radii)*cosd(45),'-','FontSize',RecLabelFS,...
    'FontSize',RecLabelFS,'HorizontalAlignment','Center','VerticalAlignment','Middle')
text(0.5 - (mean(DP6.radii) + DP6.radii(2))*cosd(45),0.5 + (mean(DP6.radii) + DP6.radii(2))*cosd(45),'No Change',...
    'FontSize',RecLabelFS)

% Plot the stimulation locations
for i = 1:3
    scatter(DP6.stimLocs(i,1),DP6.stimLocs(i,2),StimDotSz,PColor{i},'filled')
    text(DP6.stimLocs(i,1) + 0.02,DP6.stimLocs(i,2),['Stim ',num2str(i)],'FontSize',StimDotFS,'Color',PColor{i})
end

xlim([0,1])
ylim([0,1])
set(gca,'FontSize',UnitFS)
xlabel('x position','FontSize',AxLabelFS)
ylabel('y position','FontSize',AxLabelFS)
title('Example Receptive Field','FontSize',TitleFS)



% Plot the example spike responses

f = subplot('Position',[LeftCoord1(2),BottomCoord(1),width2,height1]);
hold on

fill([DP6.timeLocs(2),DP6.timeLocs(3),DP6.timeLocs(3),DP6.timeLocs(2)],[0,0,1,1],[0.8,0.8,0.8],'EdgeColor','None')

% Set the parameters for the vertical space
yBot = 0.05;
yTop = 0.15;
yBuff = 0.1;
ySpk = (1 - yBot - yTop - yBuff*2)/3;

yVals = [1 - yTop,1 - yTop - ySpk,1 - yTop - ySpk - yBuff,1 - yTop - 2*ySpk - yBuff,1 - yTop - 2*ySpk - 2*yBuff, 1 - yTop - 3*ySpk - 2*yBuff];

for iSpike = 1:length(DP6.spikes{1})
    line([DP6.spikes{1}(iSpike),DP6.spikes{1}(iSpike)],[yVals(1),yVals(2)],'Color',PColor{1},'LineWidth',SpikeLW)
end
for iSpike = 1:length(DP6.spikes{2})
    line([DP6.spikes{2}(iSpike),DP6.spikes{2}(iSpike)],[yVals(3),yVals(4)],'Color',PColor{2},'LineWidth',SpikeLW)
end
for iSpike = 1:length(DP6.spikes{3})
    line([DP6.spikes{3}(iSpike),DP6.spikes{3}(iSpike)],[yVals(5),yVals(6)],'Color',PColor{3},'LineWidth',SpikeLW)
end

text(mean(DP6.timeLocs(2:3)),mean([1,1 - yTop]),'Stimulation On','FontSize',StimLabelFS,'HorizontalAlignment','Center')
text(mean(DP6.timeLocs(1:2)),yVals(1) + 0.01,'Stim 1','Color',PColor{1},'FontSize',StimLabelFS,'VerticalAlignment','Bottom','HorizontalAlignment','Center')
text(mean(DP6.timeLocs(1:2)),yVals(3) + 0.01,'Stim 2','Color',PColor{2},'FontSize',StimLabelFS,'VerticalAlignment','Bottom','HorizontalAlignment','Center')
text(mean(DP6.timeLocs(1:2)),yVals(5) + 0.01,'Stim 3','Color',PColor{3},'FontSize',StimLabelFS,'VerticalAlignment','Bottom','HorizontalAlignment','Center')

xlim([DP6.timeLocs(1),DP6.timeLocs(4)])
ylim([0,1])

set(gca,'YTickLabel',{})
set(gca,'YTick',[])
set(gca,'FontSize',UnitFS)
set(gca,'Layer','top')
xlabel('Time (ms)','FontSize',AxLabelFS)
title('Example Spike Responses','FontSize',TitleFS)




% Plot the stimulus position encoding

f = subplot('Position',[LeftCoord1(3),BottomCoord(1),width2,height1]);
hold on

yMax = max(DP7.errorBarVals(:));
yMin = min(DP7.errorBarVals(:));
dy = yMax - yMin;

fill([DP6.timeLocs(2),DP6.timeLocs(3),DP6.timeLocs(3),DP6.timeLocs(2)],[yMin - 0.1*dy,yMin - 0.1*dy,yMax + 0.3*dy,yMax + 0.3*dy],[0.8,0.8,0.8],'EdgeColor','None')
text(mean(DP6.timeLocs(2:3)),mean([yMax + 0.1*dy,yMax + 0.3*dy]),'Stimulation On','FontSize',StimLabelFS,'HorizontalAlignment','Center')

for iT = 1:size(DP7.errorBarVals,1)
    scatter(DP7.time(iT),DP7.errorBarVals(iT,1),MIDotSz,'k','filled')
    line([DP7.time(iT),DP7.time(iT)],[DP7.errorBarVals(iT,2),DP7.errorBarVals(iT,3)],'LineWidth',ErrBarLW,'Color','k')
end

ylim([yMin - 0.1*yMin,yMax + 0.3*dy])
xlim(DP7.xLims)


set(gca,'FontSize',UnitFS)
set(gca,'Layer','top')
xlabel('Time (ms)','FontSize',AxLabelFS)
ylabel('Mutual Information (bits)','FontSize',AxLabelFS)
title('Stimulus Location Encoding','FontSize',TitleFS)




% Plot the PID values as a function of distance between place fields

yMax = max(max(max(DP15.PIDCurves(:,:,1:3))));
yMin = min(min(min(DP15.PIDCurves(:,:,1:3))));
dy = yMax - yMin;

for iMeas = 1:4
    
    f = subplot('Position',[LeftCoord2(iMeas),BottomCoord(2),width3,height2]);
    hold on
    
    scatter(DP15.dist,DP15.PIDCurves(:,iMeas,1),PIDDotSz,'k','filled')
    for iDist = 1:length(DP14.dist)
        line([DP15.dist(iDist),DP15.dist(iDist)],[DP15.PIDCurves(iDist,iMeas,2),DP15.PIDCurves(iDist,iMeas,3)],'Color','k','LineWidth',PIDLW)
    end
    
    xlim([0,DP15.dist(end)])
    ylim([yMin - 0.1*dy,yMax + 0.1*dy])
    
    set(gca,'FontSize',UnitFS)
    set(gca,'Layer','Top')
    
    xlabel('Receptive Field Separation','FontSize',AxLabelFS)
    ylabel(PIDNamesShort{iMeas},'FontSize',AxLabelFS)
    title(PIDNamesLong{iMeas},'FontSize',TitleFS)
    
end



% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo5CS')
cd(SimDir)






%% Make the Place Cell Figure

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.4;

% Set the top and bottom margins in inches
tmargin = 0.25;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace1 = 0.1; % space between heat maps on left and color bars on right (or between two heat maps)
hspace2 = 0.7; % space between color bar and heat map (with axis labels)
hspace3 = 0.4; % space between diagram and heat map
hspace4 = 0.4; % space between color bar and box plot
hspace5 = 0.5; % space between distance dependent PID plots

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set the width of the colorbar in inches
widthCB = 0.2;

% Set the width of the box plot
widthBP = 0.5;

% Set Paper size
papersize = [7.4 6];

% calculate the width and height based on dimensions above in inches
width1 = (papersize(1) - lmargin - hspace3 - rmargin - 2*hspace1 - hspace2 - 2*widthCB - widthBP - hspace4)/3;
width2 = (papersize(1) - lmargin - (rmargin - 0.3) - 3*hspace5)/4;
height1 = width1;
height2 = width2;

LeftCoord1 = [lmargin,...
    lmargin + width1 + hspace3,...
    lmargin + width1 + hspace3 + width1 + hspace1,...
    lmargin + width1 + hspace3 + width1 + hspace1 + widthCB + hspace2,...
    lmargin + width1 + hspace3 + width1 + hspace1 + widthCB + hspace2 + width1 + hspace1,...
    lmargin + width1 + hspace3 + width1 + hspace1 + widthCB + hspace2 + width1 + hspace1 + widthCB + hspace4];
% LeftCoord1 = [lmargin,2*lmargin + width1,2*lmargin + 2*width1 + hspace1,2*lmargin + 2*width1 + hspace1 + widthCB + hspace2,2*lmargin + 3*width1 + 2*hspace1 + widthCB + hspace2];
LeftCoord2 = (0:(4 - 1))*(width2 + hspace5) + lmargin;
%LeftCoord2 = (0:(5 - 1))*(width2 + hspace1) + lmargin;
BottomCoord = [papersize(2) - tmargin - height1,papersize(2) - tmargin - height1 - vspace - height2];

% Convert to fraction of page sizes
LeftCoord1 = LeftCoord1/papersize(1);
LeftCoord2 = LeftCoord2/papersize(1);
BottomCoord = BottomCoord/papersize(2);
width1 = width1/papersize(1);
width2 = width2/papersize(1);
widthCB = widthCB/papersize(1);
widthBP = widthBP/papersize(1);
height1 = height1/papersize(2);
height2 = height2/papersize(2);


% Set Colors
PColor = {[0,0,1],[1,0,0],[0,210/255,50/255],[150/255,0,200/255],[238/255,130/255,238/255],[200/255,0,255/255],[0/255,255/255,185/255],[238/255,130/255,238/255]};
PColorGS = cell([5,1]);
for i = 1:5
    PColorGS{i} = [0,0,0] + 0.15*(i - 1);
end

% Set the background gray color
BGColor = 0.85*ones([1,3]);

% Set fontsizes
LegFS = 6;
TitleFS = 8;
AxLabelFS = 7;
UnitFS = 6;
StimLabelFS = 8;
AngleNameFS = 8;
StimDotFS = 8;
RecLabelFS = 8;
CBFS = 6;

% Set the linewidths
VoltageLW = 0.6;
ThickLW = 5;
ThinLW = 1;
StimLineLW = 0.8;
DelayLW = 2;
AxLW = 0.4;
SpikeLW = 0.4;
FRLW = 1;
DirLW = 5;
MIRespLW = 0.8;
CircleLW = 1.2;
ErrBarLW = 0.8;
PlaceDotLW = 1.2;
PIDLW = 1;

% Set dot sizes
StimDotSz = 20;
MIDotSz = 15;
BPDotSz = 20;
PlaceDotSz = 20;

% Set the significance dot offset and size
sigDotRat = 0.08;
sigDotSz = 15;

% Make PID names
PIDNamesLong = {'Redundancy','Unique Center','Unique Variable','Synergy'};
PIDNamesShort = {'Red (bits)','Uniq A (bits)','Uniq B (bits)','Syn (bits)'};

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');



% Make the diagram

f = subplot('Position',[LeftCoord1(1),BottomCoord(1),width1,height1]);
hold on

xlim([0,1])
ylim([0,1])
ticks = 0:0.2:1;
ticklabels = cell([1,length(ticks)]);
for i = 1:length(ticks)
    ticklabels{i} = num2str(ticks(i));
end

set(gca,'FontSize',UnitFS)
set(gca,'YTick',ticks)
set(gca,'YTickLabel',ticklabels)
set(gca,'XTick',ticks)
set(gca,'XTickLabel',ticklabels)
xlabel('x position','FontSize',AxLabelFS)
ylabel('y position','FontSize',AxLabelFS)

title('Place Cells in an Exploring Animal','FontSize',TitleFS)



% Make the linger time heat map

f = subplot('Position',[LeftCoord1(2),BottomCoord(1),width1,height1]);
hold on

cMin = min(DP9.anHM(:));
cMax = max(DP9.anHM(:));

nBins = size(DP9.anHM,1);

imagesc(DP9.anHM',[cMin,cMax])

xlim([0.5,nBins + 0.5])
ylim([0.5,nBins + 0.5])

m = nBins;
b = 0.5;
ticks = 0:0.2:1;
ticklabels = cell([1,length(ticks)]);
for i = 1:length(ticks)
    ticklabels{i} = num2str(ticks(i));
end

set(gca,'FontSize',UnitFS)
set(gca,'Layer','top')
set(gca,'YTick',m*ticks + b)
set(gca,'YTickLabel',ticklabels)
set(gca,'XTick',m*ticks + b)
set(gca,'XTickLabel',ticklabels)
xlabel('x position','FontSize',AxLabelFS)
ylabel('y position','FontSize',AxLabelFS)

set(gca,'YDir','normal')

title('Example Animal Linger Time','FontSize',TitleFS)

% Make the colorbar
hbar = colorbar;
set(hbar,'location','manual',...
    'position',[LeftCoord1(3),BottomCoord(1),widthCB,height1],...
    'fontsize',CBFS,...
    'ylim',[cMin,cMax]);
ylabel(hbar,'Time (s)','FontSize',AxLabelFS)





% Make the firing rate time heat map

f = subplot('Position',[LeftCoord1(4),BottomCoord(1),width1,height1]);
hold on

cMin = min(DP9.spkHM(:));
cMax = max(DP9.spkHM(:));

nBins = size(DP9.spkHM,1);

imagesc(DP9.spkHM',[cMin,cMax])

xlim([0.5,nBins + 0.5])
ylim([0.5,nBins + 0.5])

m = nBins;
b = 0.5;
ticks = 0:0.2:1;
ticklabels = cell([1,length(ticks)]);
for i = 1:length(ticks)
    ticklabels{i} = num2str(ticks(i));
end

% Add the place dot
scatter(m*0.5 + b,m*0.5 + b,PlaceDotSz,'w','LineWidth',PlaceDotLW)

set(gca,'FontSize',UnitFS)
set(gca,'Layer','top')
set(gca,'YTick',m*ticks + b)
set(gca,'YTickLabel',ticklabels)
set(gca,'XTick',m*ticks + b)
set(gca,'XTickLabel',ticklabels)
xlabel('x position','FontSize',AxLabelFS)
ylabel('y position','FontSize',AxLabelFS)

set(gca,'YDir','normal')

title('Example Neuron Firing Rate','FontSize',TitleFS)

% Make the colorbar
hbar = colorbar;
set(hbar,'location','manual',...
    'position',[LeftCoord1(5),BottomCoord(1),widthCB,height1],...
    'fontsize',CBFS,...
    'ylim',[cMin,cMax]);
ylabel(hbar,'Firing Rate (Hz)','FontSize',AxLabelFS)



% Make the box plot

f = subplot('Position',[LeftCoord1(6),BottomCoord(1),widthBP,height1]);
hold on

yMax = max(DP10.errorBarVals(5,:));
yMin = min(DP10.errorBarVals(1,:));
dy = yMax - yMin;
if DP10.sigDots > 0
    yMax = yMax + dy*sigDotRat*(DP10.sigDots + 1);
end

line([1,1],[DP10.errorBarVals(1,1),DP10.errorBarVals(5,1)],'Color','k','LineWidth',ThinLW)
line([1,1],[DP10.errorBarVals(2,1),DP10.errorBarVals(4,1)],'Color','k','LineWidth',ThickLW)
line([2,2],[DP10.errorBarVals(1,2),DP10.errorBarVals(5,2)],'Color','k','LineWidth',ThinLW)
line([2,2],[DP10.errorBarVals(2,2),DP10.errorBarVals(4,2)],'Color','k','LineWidth',ThickLW)

if DP10.sigDots > 0
    xCoords = 1.5*ones([1,DP10.sigDots]);
    yCoords = max(DP10.errorBarVals(5,:)) + dy*sigDotRat*(1:DP10.sigDots);
    scatter(xCoords,yCoords,sigDotSz,'k','filled')
end

xlim([0,3])
ylim([yMin - 0.03*dy,yMax + 0.03*dy])

set(gca,'FontSize',UnitFS)
set(gca,'XTick',[1,2])
set(gca,'XTickLabel',[])
set(gca,'YAxisLocation','Right')
ylabel('Mutual Information (bits)','FontSize',AxLabelFS)
title('Location Encoding','FontSize',TitleFS)




% Plot the PID values as a function of distance between place fields

yMax = max(max(max(DP14.PIDCurves(:,:,1:3))));
yMin = min(min(min(DP14.PIDCurves(:,:,1:3))));
dy = yMax - yMin;

for iMeas = 1:4
    
    f = subplot('Position',[LeftCoord2(iMeas),BottomCoord(2),width2,height2]);
    hold on
    
    scatter(DP14.dist,DP14.PIDCurves(:,iMeas,1),PIDDotSz,'k','filled')
    for iDist = 1:length(DP14.dist)
        line([DP14.dist(iDist),DP14.dist(iDist)],[DP14.PIDCurves(iDist,iMeas,2),DP14.PIDCurves(iDist,iMeas,3)],'Color','k','LineWidth',PIDLW)
    end
    
    xlim([0,DP14.dist(end)])
    ylim([yMin - 0.1*dy,yMax + 0.1*dy])
    
    set(gca,'FontSize',UnitFS)
    set(gca,'Layer','Top')
    
    xlabel('Place Field Separation','FontSize',AxLabelFS)
    ylabel(PIDNamesShort{iMeas},'FontSize',AxLabelFS)
    title(PIDNamesLong{iMeas},'FontSize',TitleFS)
    
end



% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo5PC')
cd(SimDir)








%% Make the Sensory Habituation Figure

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.2;

% Set the top and bottom margins in inches
tmargin = 0.25;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace1 = 0.4;

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set Paper size
papersize = [7.4 4];

% calculate the width and height based on dimensions above in inches
width1 = (papersize(1) - lmargin - rmargin - hspace1)/2;
width2 = (papersize(1) - lmargin - rmargin - 3*hspace1)/4;
height = (papersize(2) - tmargin - bmargin - vspace)/2;

LeftCoord1 = (0:(2 - 1))*(width1 + hspace1) + lmargin;
LeftCoord2 = (0:(4 - 1))*(width2 + hspace1) + lmargin;
BottomCoord = [papersize(2) - tmargin - height,bmargin];

% Convert to fraction of page sizes
LeftCoord1 = LeftCoord1/papersize(1);
LeftCoord2 = LeftCoord2/papersize(1);
BottomCoord = BottomCoord/papersize(2);
width1 = width1/papersize(1);
width2 = width2/papersize(1);
height = height/papersize(2);

% Set Colors
PColor = {[1,0,0],[0,0,1],[0,210/255,50/255],[150/255,0,200/255],[238/255,130/255,238/255],[200/255,0,255/255],[0/255,255/255,185/255],[238/255,130/255,238/255]};
PColorGS = cell([5,1]);
for i = 1:5
    PColorGS{i} = [0,0,0] + 0.15*(i - 1);
end

% Set the background gray color
BGColor = 0.85*ones([1,3]);

% Set fontsizes
LegFS = 6;
TitleFS = 8;
AxLabelFS = 7;
UnitFS = 6;
StimLabelFS = 6;
AngleNameFS = 8;
StimDotFS = 6;
RecLabelFS = 8;
CBFS = 6;

% Set the linewidths
VoltageLW = 0.6;
ThickLW = 5;
ThinLW = 1;
StimLineLW = 0.8;
DelayLW = 2;
AxLW = 0.4;
SpikeLW = 0.4;
FRLW = 1;
DirLW = 5;
MIRespLW = 0.8;
CircleLW = 1.2;
ErrBarLW = 0.8;
PlaceDotLW = 1.2;

% Set dot sizes
StimDotSz = 20;
MIDotSz = 15;
BPDotSz = 20;
PlaceDotSz = 20;

% Make PID names
PIDNamesLong = {'Redundancy','Unique Center','Unique Variable','Synergy'};
PIDNamesShort = {'Red (bits)','Uniq A (bits)','Uniq B (bits)','Syn (bits)'};

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');



% Make the diagram

f = subplot('Position',[LeftCoord1(1),BottomCoord(1),width1,height]);
hold on


title('Stimulus Response Habituation','FontSize',TitleFS)
axis off



% Make the example spike raster plot

f = subplot('Position',[LeftCoord1(2),BottomCoord(1),width1,height]);
hold on

fill([DP12.timeLocs(4),DP12.timeLocs(5),DP12.timeLocs(5),DP12.timeLocs(4)],[0,0,1,1],[0.8,0.8,0.8],'EdgeColor','None')

% Set the parameters for the vertical space
yBot = 0.05;
yTop = 0.15;
yBuff1 = 0.1;
yBuff2 = 0.06;
ySpk = (1 - yBot - yTop - 2*yBuff1 - 3*yBuff2)/6;

yVals = zeros([3,2,2]);
yVals(3,2,2) = yBot;
yVals(3,2,1) = yBot + ySpk;
yVals(3,1,2) = yBot + ySpk + yBuff2;
yVals(3,1,1) = yBot + 2*ySpk + yBuff2;
yVals(2,:,:) = yVals(3,:,:) + 2*ySpk + yBuff1 + yBuff2;
yVals(1,:,:) = yVals(2,:,:) + 2*ySpk + yBuff1 + yBuff2;

for iSet = 1:3
    for iNeuron = 1:2
        for iSpike = 1:length(DP12.spikes{iSet,iNeuron})
            line([DP12.spikes{iSet,iNeuron}(iSpike),DP12.spikes{iSet,iNeuron}(iSpike)],[yVals(iSet,iNeuron,1),yVals(iSet,iNeuron,2)],'Color',PColor{iNeuron},'LineWidth',SpikeLW)
        end
    end
end

for iSet = 1:3
    text(mean(DP12.timeLocs(1:2)),yVals(iSet,1,1),'S Neuron','Color',PColor{1},'FontSize',StimLabelFS,'VerticalAlignment','Bottom','HorizontalAlignment','Left')
    text(mean(DP12.timeLocs(1:2)),yVals(iSet,2,1),'M Neuron','Color',PColor{2},'FontSize',StimLabelFS,'VerticalAlignment','Bottom','HorizontalAlignment','Left')
end
text(mean(DP12.timeLocs(4:5)),mean([1,1 - yTop]),'Stimulation On','FontSize',StimLabelFS,'HorizontalAlignment','Center')

xlim([DP12.timeLocs(1),DP12.timeLocs(end)])
ylim([0,1])
set(gca,'FontSize',UnitFS)
set(gca,'Layer','Top')

xlabel('Time (ms)','FontSize',UnitFS)
set(gca,'YTick',[mean([yVals(3,1,2),yVals(3,2,1)]),mean([yVals(2,1,2),yVals(2,2,1)]),mean([yVals(1,1,2),yVals(1,2,1)])])
set(gca,'YTickLabel',{'Last Trial','Middle Trial','First Trial'})
title('Example Spike Rasters','FontSize',TitleFS)



% Plot the information results

yMax = max(DP13.errorBarVals(:));
yMin = min(DP13.errorBarVals(:));
dy = yMax - yMin;

for iMeas = 1:4
    
    f = subplot('Position',[LeftCoord2(iMeas),BottomCoord(2),width2,height]);
    hold on
    
    fill([DP12.timeLocs(2),DP12.timeLocs(3),DP12.timeLocs(3),DP12.timeLocs(2)],[yMin - 0.1*dy,yMin - 0.1*dy,yMax + 0.3*dy,yMax + 0.3*dy],[0.8,0.8,0.8],'EdgeColor','None')
    text(mean(DP12.timeLocs(2:3)),mean([yMax + 0.1*dy,yMax + 0.3*dy]),'Stimulation On/Off','FontSize',StimLabelFS,'HorizontalAlignment','Center')
    
    for iT = 1:size(DP13.errorBarVals,1)
        scatter(DP13.time(iT),DP13.errorBarVals(iT,iMeas,1),MIDotSz,'k','filled')
        line([DP13.time(iT),DP13.time(iT)],[DP13.errorBarVals(iT,iMeas,2),DP13.errorBarVals(iT,iMeas,3)],'LineWidth',ErrBarLW,'Color','k')
    end
    
    ylim([yMin - 0.1*dy,yMax + 0.3*dy])
    xlim(DP7.xLims)
    
    set(gca,'FontSize',UnitFS)
    set(gca,'Layer','top')
    xlabel('Time (ms)','FontSize',AxLabelFS)
    ylabel('Mutual Information (bits)','FontSize',AxLabelFS)
    if iMeas == 1
        title('S Stim Encoding','FontSize',TitleFS)
    elseif iMeas == 2
        title('M Stim Encoding','FontSize',TitleFS)
    elseif iMeas == 3
        title('S Trial Encoding','FontSize',TitleFS)
    elseif iMeas == 4
        title('M Trial Encoding','FontSize',TitleFS)
    end
    
end



% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo5SH')
cd(SimDir)

