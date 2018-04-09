%% neuroDemo3
%
%   neuroDemo3 reproduces the Izhikevich network simulations. Note, the 
%   scratch directory in the simulations folder will be  used to save data. 
%   These simulations will take up about 2 hours to run and it will 
%   require about 0.5 GB of space.
%
% Other m-files required: modeltype3, instinfo, inverseFNoise
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% October 2016; Last revision: 31-Oct-2016

%% Set Parameters

disp('Note that this simulation will take about 2 hours to run and it requires about 0.5 GB of space.')

% Set the number of models to run for each subtype
nModelsPerSubtype = 20;

% Set the number of subtypes (this must be 2)
nSTs = 2;

% Set the bin size in milliseconds
BinSize = 5;

% Set the time window around the start of the stimuli to analyze in
% milliseconds
MaxLead = 15;
MaxLag = 70;

% Get the scratch directory to use to store data
if ispc
    ScratchDir = [pwd,'\Scratch\ModelType3\'];
elseif isunix
    ScratchDir = [pwd,'/Scratch/ModelType3/'];
elseif ismac
    ScratchDir = [pwd,'/Scratch/ModelType3/'];
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
        
        if iSubType == 1
            
            [spkTimes,neuronID,stim,conmat,modelParams,neuronPos] = modeltype3('OneStim','conMode','Spatial', 'spatialDist', 'Line');
        
        elseif iSubType == 2
            
            [spkTimes,neuronID,stim,conmat,modelParams,neuronPos] = modeltype3('TwoStim','conMode','Spatial', 'spatialDist', 'Line');
            
        end
        
        save([ScratchDir,'ModelType3ST',num2str(iSubType),'Run',num2str(iModel),'.mat'],'spkTimes','neuronID','stim','conmat','modelParams','neuronPos')

        disp(['Finished Generating Data for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Convert the Raw Data to DataRaster Format

% Go to the directory with the analysis software
cd(AnaDir)

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType3ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])
        
        % Find the number of neurons
        nNeurons = size(spkTimes,1);
        
        % Make the Data Raster
        DataRaster = cell([2,1]);
        
        if iSubType == 1
            
            % Subtype 1 has only one type of stimulus
            
            % Format the Data
            [Temp,timeboundaries] = formattool(ones(size(spkTimes{1})),spkTimes{1},stim.times,BinSize,MaxLead,MaxLag,'count');
            DataRaster{1} = NaN([nNeurons,size(Temp,2),size(Temp,3)]);
            DataRaster{1}(1,:,:) = Temp;
            for i = 2:nNeurons
                DataRaster{1}(i,:,:) = formattool(ones(size(spkTimes{i})),spkTimes{i},stim.times,BinSize,MaxLead,MaxLag,'count');
            end
            DataRaster{2} = reshape(stim.types,[1,1,length(stim.types)]);
            
        elseif iSubType == 2
            
            % Subtype 2 has two types of stimuli
            
            % Check that the stimulation times are identical
            if ~isequal(stim.times{1},stim.times{2})
                error('stimulation times are not equal')
            end
            
            % Format the Data
            [Temp,timeboundaries] = formattool(ones(size(spkTimes{1})),spkTimes{1},stim.times{1},BinSize,MaxLead,MaxLag,'count');
            DataRaster{1} = NaN([nNeurons,size(Temp,2),size(Temp,3)]);
            DataRaster{1}(1,:,:) = Temp;
            for i = 2:nNeurons
                DataRaster{1}(i,:,:) = formattool(ones(size(spkTimes{i})),spkTimes{i},stim.times{1},BinSize,MaxLead,MaxLag,'count');
            end
            DataRaster{2} = NaN([2,1,size(Temp,3)]);
            for i = 1:2
                DataRaster{2}(i,1,:) = reshape(stim.types{i},[1,1,length(stim.types{i})]);
            end
            
        end
        
        % Record the binsize
        modelParams.binsize = BinSize;
        modelParams.maxlead = MaxLead;
        modelParams.maxlag = MaxLag;
        modelParams.Time = timeboundaries;
        
        % Save the Data
        save([ScratchDir,'ModelType3ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'],'DataRaster','modelParams','conmat','neuronID','stim','neuronPos')
        
        disp(['Finished Converting Data to DataRaster Format for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Convert the DataRaster Format Data to StatesRasters

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType3ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])
        
        % Make the stating method assignment array (making sure to count the
        % number of pulses that were used)
        nPulses = size(modelParams.pulseSize,2);
        MethodAssign = cell([size(DataRaster{1},1) + size(DataRaster{2},1),4]);
        for i = 1:size(DataRaster{2},1)
            MethodAssign{i,1} = 2;
            MethodAssign{i,2} = i;
            MethodAssign{i,3} = 'Nat';
            MethodAssign{i,4} = {};
        end
        for i = 1:size(DataRaster{1},1)
            MethodAssign{i + size(DataRaster{2},1),1} = 1;
            MethodAssign{i + size(DataRaster{2},1),2} = i;
            MethodAssign{i + size(DataRaster{2},1),3} = 'UniCB';
            MethodAssign{i + size(DataRaster{2},1),4} = {nPulses};
        end
        
        % Convert the data
        StatesRaster = data2states(DataRaster,MethodAssign);
        
        % Save the Data
        save([ScratchDir,'ModelType3ST',num2str(iSubType),'Run',num2str(iModel),'SR.mat'],'StatesRaster','modelParams','conmat','neuronID','stim','neuronPos')
        
        disp(['Finished Converting Data to StatesRaster Format for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Perform the Information Theory Analysis

InfoResults = cell([nSTs,nModelsPerSubtype]);
SubTypeParams = cell([nSTs,1]);
for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType3ST',num2str(iSubType),'Run',num2str(iModel),'SR.mat'])
        
        % Save the subtype parameters. We assume all runs of a subtype used
        % the same parameters.
        if iModel == 1
            SubTypeParams{iSubType} = modelParams;
        end
        
        % Find the total number of time bins
        nT = size(StatesRaster{1},2);
        
        % Find the total number of neurons
        nNeurons = size(StatesRaster{1},1);
        
        % Do the task
        if iSubType == 1
            
            % Subtype 1 data have only one type of stimulus
            MIResults = NaN([nNeurons,nT]);
            if iModel == 1 % See notes near TE calculation below
                TEResults = NaN([nNeurons,nNeurons,nT]);
            else
                TEResults = NaN;
            end
            for iNeuron = 1:nNeurons
                
                % First, calculate the mutual information between the 
                % neuron activity and the stimulation state
                Method = 'PairMI';
                for iT = 1:nT
                    VariableIDs = {1,iNeuron,iT;2,1,1};
                    MIResults(iNeuron,iT) = instinfo(StatesRaster, Method, VariableIDs);
                end
            
                % Second, calculate the TE from all other neurons (only
                % need this for one example)
                if iModel == 1
                    Method = 'TE';
                    for jNeuron = setdiff(1:nNeurons,iNeuron)
                        for iT = 2:nT
                            VariableIDs = {1,iNeuron,iT;1,iNeuron,iT - 1;1,jNeuron,iT - 1};
                            TEResults(jNeuron,iNeuron,iT) = instinfo(StatesRaster, Method, VariableIDs);
                        end
                    end
                    % The previous 7 lines can be commented out and 
                    % replaced with the following line to speed up the 
                    % analysis for testing purposes
%                     TEResults(:,iNeuron,:) = rand([nNeurons,1,nT]);
                end
            end
            
            % Put the results in the right variable
            InfoResults{iSubType,iModel} = {MIResults,TEResults};
            
            
        elseif iSubType == 2
            
            % Subtype 2 data have two stimuli
            
            PIDVals = NaN([nNeurons,4,nT]);
            Method = '2PID';
            for iNeuron = 1:nNeurons
                for iT = 1:nT
                    VariableIDs = {1,iNeuron,iT;2,1,1;2,2,1};
                    PIDVals(iNeuron,:,iT) = instinfo(StatesRaster, Method, VariableIDs);
                end
            end
            
            % Put the results in the right variable
            InfoResults{iSubType,iModel} = {PIDVals};
            
        end
        
        if (iSubType == 1) && (iModel == 1)
            disp(['Finished Analyzing for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
            disp('    Note: This run took much longer to analyze than the others will due to example TE calculations.')
        else
            disp(['Finished Analyzing for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        end
        
    end
end

% Save the Results
save([ScratchDir,'ModelType3InfoResults.mat'],'InfoResults','SubTypeParams')

%% Gather Data for Plotting

% Get an Example Network (Subtype 1) for Diagram

DP1 = struct;

% Load the first run for subtype 1
load([ScratchDir,'ModelType3ST1Run1.mat'])

% Record the connection matrix
DP1.conmat = conmat;

% Get the neuron positions
DP1.pos = cell([3,1]);
DP1.pos{1} = neuronPos(strcmp(neuronID,'1'),:);
DP1.pos{2} = neuronPos(strcmp(neuronID,'Ex'),:);
DP1.pos{3} = neuronPos(strcmp(neuronID,'In'),:);

% Record the neuron identities
DP1.type = cell([3,1]);
DP1.type{1} = find(strcmp(neuronID,'1'));
DP1.type{2} = find(strcmp(neuronID,'Ex'));
DP1.type{3} = find(strcmp(neuronID,'In'));

% Get all the positions
DP1.neuronPos = neuronPos;



% Get an Example Spike Raster (Subtype 1)

DP2 = struct;

% Load the first run for subtype 1
load([ScratchDir,'ModelType3ST1Run1.mat'])

% Make lists of the neurons
In1 = find(strcmp(neuronID,'1'));
Ex = find(strcmp(neuronID,'Ex'));
Inh = find(strcmp(neuronID,'In'));

% Reorganize the list of excitatory and inhibitory neurons based on
% distance from the line x = 0.5
distance = abs(neuronPos(Ex,1) - 0.5);
[waste,I] = sort(distance,1,'ascend');
Ex = Ex(I);
distance = abs(neuronPos(Inh,1) - 0.5);
[waste,I] = sort(distance,1,'ascend');
Inh = Inh(I);

% Make a relabelled list of the neurons for plotting
In1Re = 1:length(In1);
ExRe = (1:length(Ex)) + length(In1);
InhRe = (1:length(Inh)) + length(In1) + length(Ex);
NewList = zeros([1,1000]);
NewList(In1) = In1Re;
NewList(Ex) = ExRe;
NewList(Inh) = InhRe;

% Calculate time bounds based on the stimulations and parameters
stimTime = stim.times(find(stim.types == 2,1,'first'));
StartTime = stimTime - modelParams.interstimT - 1;
EndTime = stimTime + 2*modelParams.interstimT + modelParams.intrastimT - 1;

% Make a list of the spikes
for i = 1:length(spkTimes)
    spkTimes{i}(spkTimes{i} > EndTime) = [];
    spkTimes{i}(spkTimes{i} < StartTime) = [];
end
nSpikes = 0;
for i = 1:length(spkTimes)
    nSpikes = nSpikes + length(spkTimes{i});
end
spikes = zeros([nSpikes,3]);
iSpike = 1;
for i = 1:length(spkTimes)
    spikes(iSpike:(iSpike + length(spkTimes{i}) - 1),1) = NewList(i);
    spikes(iSpike:(iSpike + length(spkTimes{i}) - 1),2) = spkTimes{i} - stimTime;
    spikes(iSpike:(iSpike + length(spkTimes{i}) - 1),3) = i;
    iSpike = iSpike + length(spkTimes{i});
end

DP2.spikes = cell([3,1]);
DP2.spikes{1} = spikes(ismember(spikes(:,3),In1),:);
DP2.spikes{2} = spikes(ismember(spikes(:,3),Ex),:);
DP2.spikes{3} = spikes(ismember(spikes(:,3),Inh),:);

DP2.StartTime = StartTime - stimTime + 1;
DP2.EndTime = EndTime - stimTime + 1;



% Get MI Heat Maps

% Set the number of spatial bins
nBins = 20;

DP3 = struct;

% Record the number of spatial bins
DP3.nBins = nBins;

% Get the number of time steps and the time windows for each bin (should be
% the same for all models)
load([ScratchDir,'ModelType3ST1Run1SR.mat'])
nT = size(StatesRaster{1},2);
DP3.Time = modelParams.Time;

% Make the spatial edges for the spatial bins
spatialEdges = linspace(0,1,nBins + 1);

% Calculate the average MI in each spatial bin
DP3.MIHeatMap = zeros([nBins,nBins,nT]);
Counts = zeros([nBins,nBins,nT]);

iSubType = 1;
for iModel = 1:nModelsPerSubtype
    load([ScratchDir,'ModelType3ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])
    for xBin = 1:nBins
        for yBin = 1:nBins
            Group1 = find((neuronPos(:,1) >= spatialEdges(xBin)) & (neuronPos(:,1) < spatialEdges(xBin + 1)));
            Group2 = find((neuronPos(:,2) >= spatialEdges(yBin)) & (neuronPos(:,2) < spatialEdges(yBin + 1)));
            FoundNeurons = intersect(Group1,Group2);
            for iT = 1:nT
                DP3.MIHeatMap(xBin,yBin,iT) = DP3.MIHeatMap(xBin,yBin,iT) + sum(squeeze(InfoResults{iSubType,iModel}{1}(FoundNeurons,iT)));
                Counts(xBin,yBin,iT) = Counts(xBin,yBin,iT) + length(FoundNeurons);
            end
        end
    end
end

for xBin = 1:nBins
    for yBin = 1:nBins
        for iT = 1:nT
            DP3.MIHeatMap(xBin,yBin,iT) = DP3.MIHeatMap(xBin,yBin,iT)/Counts(xBin,yBin,iT);
        end
    end
end



% Get TE Heat Map Values

DP4 = struct;

% Get the number of time steps and the time windows for each bin (should be
% the same for all models)
load([ScratchDir,'ModelType3ST1Run1SR.mat'])
nT = size(StatesRaster{1},2);
DP4.Time = modelParams.Time;

% Make lists of the neurons
In1 = find(strcmp(neuronID,'1'));
Ex = find(strcmp(neuronID,'Ex'));
Inh = find(strcmp(neuronID,'In'));

% Reorganize the list of excitatory and inhibitory neurons based on
% distance from the line x = 0.5
distance = abs(neuronPos(Ex,1) - 0.5);
[waste,I] = sort(distance,1,'ascend');
Ex = Ex(I);
distance = abs(neuronPos(Inh,1) - 0.5);
[waste,I] = sort(distance,1,'ascend');
Inh = Inh(I);

% Make a relabelled list of the neurons for plotting
In1Re = 1:length(In1);
ExRe = (1:length(Ex)) + length(In1);
InhRe = (1:length(Inh)) + length(In1) + length(Ex);
NewList = zeros([1,1000]);
NewList(In1) = In1Re;
NewList(Ex) = ExRe;
NewList(Inh) = InhRe;

% Reorganize the TE results
DP4.TE = zeros([1000,1000,nT]);
iModel = 1; % This is why we only ran the TE analysis once
for iNeuron = 1:1000
    for jNeuron = setdiff(1:1000,iNeuron)
        DP4.TE(NewList(iNeuron),NewList(jNeuron),:) = squeeze(InfoResults{iSubType,iModel}{2}(iNeuron,jNeuron,:));
    end
end
DP4.TE(:,:,1) = 0;

% Get information about the neurons
DP4.neuron = cell([3,2]);
DP4.neuron{1,1} = In1;
DP4.neuron{2,1} = Ex;
DP4.neuron{3,1} = Inh;
DP4.neuron{1,2} = In1Re;
DP4.neuron{2,2} = ExRe;
DP4.neuron{3,2} = InhRe;

DP4.NewList = NewList;



% Get an Example Network (Subtype 2) for Diagram

DP5 = struct;

% Load the first run for subtype 2
load([ScratchDir,'ModelType3ST2Run1.mat'])

% Record the connection matrix
DP5.conmat = conmat;

% Get the neuron positions
DP5.pos = cell([4,1]);
DP5.pos{1} = neuronPos(strcmp(neuronID,'1'),:);
DP5.pos{2} = neuronPos(strcmp(neuronID,'2'),:);
DP5.pos{3} = neuronPos(strcmp(neuronID,'Ex'),:);
DP5.pos{4} = neuronPos(strcmp(neuronID,'In'),:);

% Record the neuron identities
DP5.type = cell([4,1]);
DP5.type{1} = find(strcmp(neuronID,'1'));
DP5.type{2} = find(strcmp(neuronID,'2'));
DP5.type{3} = find(strcmp(neuronID,'Ex'));
DP5.type{4} = find(strcmp(neuronID,'In'));

% Get all the positions
DP5.neuronPos = neuronPos;



% Get an Example Spike Raster (Subtype 2)

DP6 = struct;

% Load the first run for subtype 2
load([ScratchDir,'ModelType3ST2Run1.mat'])

% Make lists of the neurons
In1 = find(strcmp(neuronID,'1'));
In2 = find(strcmp(neuronID,'2'));
Ex = find(strcmp(neuronID,'Ex'));
Inh = find(strcmp(neuronID,'In'));

% Make a new list of the neurons to ease reorganization
NewList1 = zeros([800,3]);
NewList1(:,1) = 1:800;
NewList1(:,2) = neuronPos(1:800,1);
NewList1(In1,3) = 1;
NewList1(In2,3) = 2;
NewList1(Ex,3) = 3;
NewList2 = zeros([200,3]);
NewList2(:,1) = 801:1000;
NewList2(:,2) = neuronPos(801:1000,1);
NewList2(:,3) = 4;

% Sort the new lists based on x position
NewList1 = sortrows(NewList1,2);
NewList2 = sortrows(NewList2,2);

% Calculate time bounds based on the stimulations and parameters
stimTime = zeros([1,3]);
stimTime(1) = stim.times{1}(find((stim.types{1} == 1) & (stim.types{2} == 2),1,'first'));
stimTime(2) = stim.times{1}(find((stim.types{1} == 2) & (stim.types{2} == 1),1,'first'));
stimTime(3) = stim.times{1}(find((stim.types{1} == 2) & (stim.types{2} == 2),1,'first'));
StartTime = stimTime - modelParams.interstimT - 1;
EndTime = stimTime + 2*modelParams.interstimT + modelParams.intrastimT - 1;

% Make a list of the spikes
nSpikes = zeros([1,3]);
spikes = cell([3,1]);
for iStim = 1:3
    for iNeuron = 1:length(spkTimes)
        nSpikes(iStim) = nSpikes(iStim) + length(spkTimes{iNeuron}((spkTimes{iNeuron} >= StartTime(iStim)) & (spkTimes{iNeuron} <= EndTime(iStim))));
    end
    spikes{iStim} = zeros([nSpikes(iStim),3]);
    iSpike = 1;
    for iNeuron = 1:length(spkTimes)
        temp = spkTimes{iNeuron}((spkTimes{iNeuron} >= StartTime(iStim)) & (spkTimes{iNeuron} <= EndTime(iStim)));
        if iNeuron <= 800
            spikes{iStim}(iSpike:(iSpike + length(temp) - 1),1) = find(NewList1(:,1) == iNeuron);
            spikes{iStim}(iSpike:(iSpike + length(temp) - 1),2) = temp - stimTime(iStim);
            spikes{iStim}(iSpike:(iSpike + length(temp) - 1),3) = NewList1(NewList1(:,1) == iNeuron,3);
        else
            spikes{iStim}(iSpike:(iSpike + length(temp) - 1),1) = find(NewList2(:,1) == iNeuron) + 800;
            spikes{iStim}(iSpike:(iSpike + length(temp) - 1),2) = temp - stimTime(iStim);
            spikes{iStim}(iSpike:(iSpike + length(temp) - 1),3) = NewList2(NewList2(:,1) == iNeuron,3);
        end
        iSpike = iSpike + length(temp);
    end    
end

DP6.spikes = cell([4,3]);
for iStim = 1:3
    DP6.spikes{1,iStim} = spikes{iStim}(spikes{iStim}(:,3) == 1,:);
    DP6.spikes{2,iStim} = spikes{iStim}(spikes{iStim}(:,3) == 2,:);
    DP6.spikes{3,iStim} = spikes{iStim}(spikes{iStim}(:,3) == 3,:);
    DP6.spikes{4,iStim} = spikes{iStim}(spikes{iStim}(:,3) == 4,:);
end

DP6.StartTime = StartTime - stimTime + 1;
DP6.EndTime = EndTime - stimTime + 1;



% Get PID Heat Maps

% Set the number of spatial bins
nBins = 20;

DP7 = struct;

% Record the number of spatial bins
DP7.nBins = nBins;

% Get the number of time steps and the time windows for each bin (should be the same for all models)
load([ScratchDir,'ModelType3ST2Run1SR.mat'])
nT = size(StatesRaster{1},2);
DP7.Time = modelParams.Time;

% Make the spatial edges for the spatial bins
spatialEdges = linspace(0,1,nBins + 1);

% Calculate the average MI in each spatial bin
DP7.PIDHeatMap = zeros([nBins,nBins,nT,4]);
Counts = zeros([nBins,nBins,nT,4]);

iSubType = 2;
for iModel = 1:nModelsPerSubtype
    load([ScratchDir,'ModelType3ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])
    for xBin = 1:nBins
        for yBin = 1:nBins
            Group1 = find((neuronPos(:,1) >= spatialEdges(xBin)) & (neuronPos(:,1) < spatialEdges(xBin + 1)));
            Group2 = find((neuronPos(:,2) >= spatialEdges(yBin)) & (neuronPos(:,2) < spatialEdges(yBin + 1)));
            FoundNeurons = intersect(Group1,Group2);
            for iT = 1:nT
                for iPID = 1:4
                    DP7.PIDHeatMap(xBin,yBin,iT,iPID) = DP7.PIDHeatMap(xBin,yBin,iT,iPID) + sum(squeeze(InfoResults{iSubType,iModel}{1}(FoundNeurons,iPID,iT)));
                    Counts(xBin,yBin,iT,iPID) = Counts(xBin,yBin,iT,iPID) + length(FoundNeurons);
                end
            end
        end
    end
end

DP7.PIDHeatMap = DP7.PIDHeatMap./Counts;

% Save the results
save([ScratchDir,'ModelType3PlotData'],'DP1','DP2','DP3','DP4','DP5','DP6','DP7')


%% Make the MI and TE figure

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.4;

% Set the top and bottom margins in inches
tmargin = 0.3;
bmargin = 0.4;

% Set the horizontal distance between plots in inches
hspace = [0.6,0.1];

% Set the vertical distance between plots in inches
vspace = 0.55;

% Set the width of the color bar in inches
cbwidth = 0.2;

% Set Paper size
papersize = [7 8];

% Set the number of heat maps
nHM = 5;

% calculate the width and height
width1 = (papersize(1) - lmargin - rmargin - hspace(1))/2;
width2 = (papersize(1) - lmargin - rmargin - nHM*hspace(2) - cbwidth)/nHM;
height2 = width2;
height1 = width1;

LeftCoord1 = [lmargin,lmargin + width1 + hspace(1)];
LeftCoord2 = [(0:(nHM - 1))*(width2 + hspace(2)) + lmargin,nHM*(width2 + hspace(2)) + lmargin];
BottomCoord1 = papersize(2) - tmargin - height1;
BottomCoord2 = papersize(2) - tmargin - height1 - vspace - height2 - [0,vspace + height2];

% Convert to fraction of page sizes
LeftCoord1 = LeftCoord1/papersize(1);
BottomCoord1 = BottomCoord1/papersize(2);
width1 = width1/papersize(1);
height1 = height1/papersize(2);
LeftCoord2 = LeftCoord2/papersize(1);
BottomCoord2 = BottomCoord2/papersize(2);
width2 = width2/papersize(1);
height2 = height2/papersize(2);
cbwidth = cbwidth/papersize(1);


% Set Colors
PColor = {[0,0,1],[1,0,0],[0,210/255,50/255],[0.6,0.6,0.6],[238/255,130/255,238/255],[200/255,0,255/255],[0/255,255/255,185/255],[238/255,130/255,238/255]};
PColorGS = cell([5,1]);
for i = 1:5
    PColorGS{i} = [0,0,0] + 0.15*(i - 1);
end

% Set the backgroun gray color
BGColor = 0.85*ones([1,3]);

% Set fontsizes
LegFS = 6;
TitleFS = 8;
AxLabelFS = 7;
UnitFS = 6;
CBFS = 6;

% Set the linewidths
NeuronLineLW = 0.7;
ConLW = 0.5;
StimTLW = 1;

% Set dot sizes
neurDotSz = 10;
spikeDotSz = 1;

% Set the time bins to use for the heat maps (there should be nHM)
HMTBins = [4,5,6,7,8];

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');



% Make the connectivity diagram (only show connections from stim neurons

f = subplot('Position',[LeftCoord1(1),BottomCoord1(1),width1,height1]);
hold on

mirroradjust = zeros([9,2]);
imirror = 1;
for i  = [-1,0,1];
    for j = [-1,0,1];
        mirroradjust(imirror,:) = [i,j];
        imirror = imirror + 1;
    end
end
% Plot the connections from stimulated neurons
for iNeuron = DP1.type{1}'
    TargetList = find(DP1.conmat(iNeuron,:) ~= 0);
    TargetList(TargetList == iNeuron) = []; % Do not plot self connections
    mirrori = repmat(DP1.neuronPos(iNeuron,:),[9,1]) + mirroradjust;
    for jNeuron = TargetList
        mirrorj = repmat(DP1.neuronPos(jNeuron,:),[9,1]) + mirroradjust;
        distance = sqrt(sum((repmat(DP1.neuronPos(iNeuron,:),[9,1]) - mirrorj).^2,2));
        minD = find(distance == min(distance),1,'first');
        if minD ~= 5
            % We must connect to a mirror
            
            % Connect from the iNeuron in the main to the jNeuron in the
            % mirror
            line([DP1.neuronPos(iNeuron,1),mirrorj(minD,1)],[DP1.neuronPos(iNeuron,2),mirrorj(minD,2)],'Color',PColor{4},'LineWidth',ConLW)
            
            % Connect from the jNeuron in the main to the iNeuron in the
            % mirror
            minD = 10 - minD;
            line([mirrori(minD,1),DP1.neuronPos(jNeuron,1)],[mirrori(minD,2),DP1.neuronPos(jNeuron,2)],'Color',PColor{4},'LineWidth',ConLW)
        else
            % We can connect without the mirror, so we only need one line
            line([DP1.neuronPos(iNeuron,1),DP1.neuronPos(jNeuron,1)],[DP1.neuronPos(iNeuron,2),DP1.neuronPos(jNeuron,2)],'Color',PColor{4},'LineWidth',ConLW)
        end
    end
end


% Plot the neuron locations
scatter(DP1.pos{2}(:,1),DP1.pos{2}(:,2),neurDotSz,PColor{1})
scatter(DP1.pos{3}(:,1),DP1.pos{3}(:,2),neurDotSz,PColor{3})
scatter(DP1.pos{1}(:,1),DP1.pos{1}(:,2),neurDotSz,PColor{2},'filled')

% Make the legend
fill([0.03,0.15,0.15,0.03],[0.03,0.03,0.16,0.16],[1,1,1],'EdgeColor','k')
x1 = 0.05;
x2 = 0.07;
y1 = 0.135;
y2 = 0.095;
y3 = 0.055;
scatter(x1,y1,neurDotSz,PColor{2},'filled')
text(x2,y1,'Stim','FontSize',LegFS)
scatter(x1,y2,neurDotSz,PColor{1})
text(x2,y2,'Ex','FontSize',LegFS)
scatter(x1,y3,neurDotSz,PColor{3})
text(x2,y3,'Inh','FontSize',LegFS)

xlabel('X Position','FontSize',AxLabelFS)
ylabel('Y Position','FontSize',AxLabelFS)
set(gca,'FontSize',UnitFS)
set(gca,'Layer','Top')
title('Example Network','FontSize',TitleFS)

xlim([0,1])
ylim([0,1])



% Make the example spike raster

f = subplot('Position',[LeftCoord1(2),BottomCoord1(1),width1,height1]);
hold on

% Set limits
xMin = -50;
xMax = 100;
yMin = 0;
yMax = 1000;


line([0,0],[0,1000],'Color','k','LineStyle','--','LineWidth',StimTLW)

% Plot the spikes
scatter(DP2.spikes{1}(:,2),DP2.spikes{1}(:,1),spikeDotSz,'filled','MarkerFaceColor',PColor{2})
scatter(DP2.spikes{2}(:,2),DP2.spikes{2}(:,1),spikeDotSz,'filled','MarkerFaceColor',PColor{1})
scatter(DP2.spikes{3}(:,2),DP2.spikes{3}(:,1),spikeDotSz,'filled','MarkerFaceColor',PColor{3})

% Make the legend
x1 = xMin + 0.83*(xMax - xMin);
x2 = xMin + 0.85*(xMax - xMin);
x3 = xMin + 0.87*(xMax - xMin);
x4 = xMin + 0.95*(xMax - xMin);
y1 = yMin + 0.03*(yMax - yMin);
y2 = yMin + 0.055*(yMax - yMin);
y3 = yMin + 0.095*(yMax - yMin);
y4 = yMin + 0.135*(yMax - yMin);
y5 = yMin + 0.16*(yMax - yMin);
fill([x1,x4,x4,x1],[y1,y1,y5,y5],[1,1,1],'EdgeColor','k')
scatter(x2,y4,neurDotSz,PColor{2},'filled')
text(x3,y4,'Stim','FontSize',LegFS)
scatter(x2,y3,neurDotSz,PColor{1},'filled')
text(x3,y3,'Ex','FontSize',LegFS)
scatter(x2,y2,neurDotSz,PColor{3},'filled')
text(x3,y2,'Inh','FontSize',LegFS)

xlim([xMin,xMax])
ylim([yMin,yMax])
set(gca,'FontSize',UnitFS)

ylabel('Neuron Number','FontSize',AxLabelFS)
xlabel('Time Relative to Stimulus (ms)','FontSize',AxLabelFS)
title('Example Spike Raster','FontSize',TitleFS)



% Make the MI heat maps


% Find the color limits
cMax = -inf;
cMin = inf;
for iPlot = 1:nHM
    cMax = max([cMax,max(max(squeeze(DP3.MIHeatMap(:,:,HMTBins(iPlot)))))]);
    cMin = min([cMin,min(min(squeeze(DP3.MIHeatMap(:,:,HMTBins(iPlot)))))]);
end
if cMin == 0
    cMin = 10^(log10(cMax) - 2);
end
    

for iPlot = 1:nHM
    f = subplot('Position',[LeftCoord2(iPlot),BottomCoord2(1),width2,height2]);
    hold on
    
    imagesc(log10(squeeze(DP3.MIHeatMap(:,:,HMTBins(iPlot)))'),log10([cMin,cMax]))
    
    line([(DP3.nBins/2) + 0.5,(DP3.nBins/2) + 0.5],[0.5,DP3.nBins + 0.5],'Color','k','LineWidth',NeuronLineLW)
    line([(DP3.nBins/2) + 0.5,(DP3.nBins/2) + 0.5],[0.5,DP3.nBins + 0.5],'Color','w','LineStyle','--','LineWidth',NeuronLineLW)
    
    xlim([0.5,DP3.nBins + 0.5])
    ylim([0.5,DP3.nBins + 0.5])
    
    set(gca,'FontSize',UnitFS)
    if iPlot > 1
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
    end
    if iPlot == 1
        set(gca,'XTick',[0.5,(DP3.nBins/2) + 0.5,DP3.nBins + 0.5])
        set(gca,'XTickLabel',{'0','0.5','1'})
        set(gca,'YTick',[0.5,(DP3.nBins/2) + 0.5,DP3.nBins + 0.5])
        set(gca,'YTickLabel',{'0','0.5','1'})
        xlabel('X Position','FontSize',AxLabelFS)
        ylabel('Y Position','FontSize',AxLabelFS)
    end
    if iPlot == round(nHM/2)
        title(['Average Stimulus Encoding (All Models)',char(10),num2str(DP3.Time(HMTBins(iPlot),1)),' to ',num2str(DP3.Time(HMTBins(iPlot),2)),' ms'],'FontSize',TitleFS)
    else
        title([num2str(DP3.Time(HMTBins(iPlot),1)),' to ',num2str(DP3.Time(HMTBins(iPlot),2)),' ms'],'FontSize',TitleFS)
    end
end

% Add the colorbar
yticks = ceil(log10(cMin)):floor(log10(cMax));
yticklabels = cell([1,length(yticks)]);
for iTick = 1:length(yticks)
    yticklabels{iTick} = num2str(yticks(iTick));
end
hbar = colorbar;
set(hbar,'location','manual',...
    'position',[LeftCoord2(nHM + 1),BottomCoord2(1),cbwidth,height2],...
    'fontsize',CBFS,...
    'YTick',yticks,...
    'YTickLabel',yticklabels,...
    'ylim',log10([cMin,cMax]));
ylabel(hbar,'log_{10} MI (bits)','FontSize',AxLabelFS)



% Make the TE heat maps

% Set the minimum value reset for plotting
minReset = 10^-5;
DP4.TE(DP4.TE < minReset) = minReset;

% Find the color limits
cMax = -inf;
cMin = inf;
for iPlot = 1:nHM
    cMax = max([cMax,max(max(squeeze(DP4.TE(:,:,HMTBins(iPlot)))))]);
    cMin = min([cMin,min(min(squeeze(DP4.TE(:,:,HMTBins(iPlot)))))]);
end

for iPlot = 1:nHM
    f = subplot('Position',[LeftCoord2(iPlot),BottomCoord2(2),width2,height2]);
    hold on
    
    imagesc(log10(flipud(squeeze(DP4.TE(:,:,HMTBins(iPlot))))),log10([cMin,cMax]))
    
    line([DP4.neuron{1,2}(end) + 0.5,DP4.neuron{1,2}(end) + 0.5],[0.5,1000 + 0.5],'Color','k','LineWidth',NeuronLineLW)
    line([DP4.neuron{2,2}(end) + 0.5,DP4.neuron{2,2}(end) + 0.5],[0.5,1000 + 0.5],'Color','k','LineWidth',NeuronLineLW)
    line([0.5,1000 + 0.5],[1000 - DP4.neuron{1,2}(end) + 0.5,1000 - DP4.neuron{1,2}(end) + 0.5],'Color','k','LineWidth',NeuronLineLW)
    line([0.5,1000 + 0.5],[1000 - DP4.neuron{2,2}(end) + 0.5,1000 - DP4.neuron{2,2}(end) + 0.5],'Color','k','LineWidth',NeuronLineLW)
    line([DP4.neuron{1,2}(end) + 0.5,DP4.neuron{1,2}(end) + 0.5],[0.5,1000 + 0.5],'Color','w','LineStyle','--','LineWidth',NeuronLineLW)
    line([DP4.neuron{2,2}(end) + 0.5,DP4.neuron{2,2}(end) + 0.5],[0.5,1000 + 0.5],'Color','w','LineStyle','--','LineWidth',NeuronLineLW)
    line([0.5,1000 + 0.5],[1000 - DP4.neuron{1,2}(end) + 0.5,1000 - DP4.neuron{1,2}(end) + 0.5],'Color','w','LineStyle','--','LineWidth',NeuronLineLW)
    line([0.5,1000 + 0.5],[1000 - DP4.neuron{2,2}(end) + 0.5,1000 - DP4.neuron{2,2}(end) + 0.5],'Color','w','LineStyle','--','LineWidth',NeuronLineLW)
    
    xlim([0.5,1000 + 0.5])
    ylim([0.5,1000 + 0.5])
    
    set(gca,'FontSize',UnitFS)
    if iPlot > 1
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
    end
    if iPlot == 1
        set(gca,'XTick',[mean([DP4.neuron{1,2}(1) - 0.5,DP4.neuron{1,2}(end) + 0.5]),...
            mean([DP4.neuron{2,2}(1) - 0.5,DP4.neuron{2,2}(end) + 0.5]),...
            mean([DP4.neuron{3,2}(1) - 0.5,DP4.neuron{3,2}(end) + 0.5])])
        set(gca,'XTickLabel',{'Stim','Ex','Inh'})
        set(gca,'YTick',fliplr([mean([1001 - DP4.neuron{1,2}(1) + 0.5,1001 - DP4.neuron{1,2}(end) - 0.5]),...
            mean([1001 - DP4.neuron{2,2}(1) + 0.5,1001 - DP4.neuron{2,2}(end) - 0.5]),...
            mean([1001 - DP4.neuron{3,2}(1) + 0.5,1001 - DP4.neuron{3,2}(end) - 0.5])]))
        set(gca,'YTickLabel',fliplr({'Stim','Ex','Inh'}))
        xlabel('Receiver','FontSize',AxLabelFS)
        ylabel('Sender','FontSize',AxLabelFS)
    end
    if iPlot == round(nHM/2)
        title(['Example Transfer Entropy',char(10),num2str(DP3.Time(HMTBins(iPlot),1)),' to ',num2str(DP3.Time(HMTBins(iPlot),2)),' ms'],'FontSize',TitleFS)
    else
        title([num2str(DP3.Time(HMTBins(iPlot),1)),' to ',num2str(DP3.Time(HMTBins(iPlot),2)),' ms'],'FontSize',TitleFS)
    end
end

% Add the colorbar
yticks = ceil(log10(cMin)):floor(log10(cMax));
yticklabels = cell([1,length(yticks)]);
for iTick = 1:length(yticks)
    yticklabels{iTick} = num2str(yticks(iTick));
end
hbar = colorbar;
set(hbar,'location','manual',...
    'position',[LeftCoord2(nHM + 1),BottomCoord2(2),cbwidth,height2],...
    'fontsize',CBFS,...
    'YTick',yticks,...
    'YTickLabel',yticklabels,...
    'ylim',log10([cMin,cMax]));
ylabel(hbar,'log_{10} TE (bits)','FontSize',AxLabelFS)


% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo3Panel1')
cd(SimDir)


%% Make the PID figure

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.45;

% Set the top and bottom margins in inches
tmargin = 0.3;
bmargin = 0.4;

% Set the horizontal distance between plots in inches
hspace = [0.59,0.15,0.12];

% Set the vertical distance between plots in inches
vspace = [0.6,0.35];

% Set the width of the color bar in inches
cbwidth = 0.2;

% Set Paper size
papersize = [7 8];

% Set the number of heat maps
nHM = 5;

% Set the width of the spike raster examples in inches
width3 = 1.225;

% calculate the width and height
width1 = papersize(1) - lmargin - rmargin - hspace(1) - 3*width3 - 2*hspace(2);
width2 = (papersize(1) - lmargin - rmargin - nHM*hspace(3) - cbwidth)/nHM;
height2 = width2;
height1 = width1;

LeftCoord1 = [lmargin,lmargin + width1 + hspace(1),lmargin + width1 + hspace(1) + width3 + hspace(2),lmargin + width1 + hspace(1) + 2*width3 + 2*hspace(2)];
LeftCoord2 = [(0:(nHM - 1))*(width2 + hspace(3)) + lmargin,nHM*(width2 + hspace(3)) + lmargin];
BottomCoord1 = papersize(2) - tmargin - height1;
BottomCoord2 = papersize(2) - tmargin - height1 - vspace(1) - height2 - [0,vspace(2) + height2,2*(vspace(2) + height2),3*(vspace(2) + height2)];

% Convert to fraction of page sizes
LeftCoord1 = LeftCoord1/papersize(1);
BottomCoord1 = BottomCoord1/papersize(2);
width1 = width1/papersize(1);
height1 = height1/papersize(2);
LeftCoord2 = LeftCoord2/papersize(1);
BottomCoord2 = BottomCoord2/papersize(2);
width2 = width2/papersize(1);
width3 = width3/papersize(1);
height2 = height2/papersize(2);
cbwidth = cbwidth/papersize(1);


% Set Some Colors
PColor = {[0,0,1],[1,0,0],[0,210/255,50/255],[0.6,0.6,0.6],[255/255,127/255,0],[200/255,0,255/255],[0/255,255/255,185/255],[238/255,130/255,238/255]};
PColorGS = cell([5,1]);
for i = 1:5
    PColorGS{i} = [0,0,0] + 0.15*(i - 1);
end

% Set the backgroun gray color
BGColor = 0.85*ones([1,3]);

% Set fontsizes
LegFS = 6;
TitleFS = 8;
AxLabelFS = 7;
UnitFS = 6;
CBFS = 6;

% Set the linewidths
NeuronLineLW = 0.7;
ConLW = 0.5;
StimTLW = 1;


% Set dot sizes
neurDotSz = 10;
spikeDotSz = 1;

% Set the PID names
PIDNamesLong = {'Redundancy','Unique A','Unique B','Synergy'};
PIDNamesShort = {'Red','Uniq A','Uniq B','Syn'};

% Set the time bins to use for the heat maps (there should be nHM)
HMTBins = [4,5,6,7,8];

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');



% Make the connectivity diagram (only show connections from stim neurons

f = subplot('Position',[LeftCoord1(1),BottomCoord1(1),width1,height1]);
hold on

mirroradjust = zeros([9,2]);
imirror = 1;
for i  = [-1,0,1];
    for j = [-1,0,1];
        mirroradjust(imirror,:) = [i,j];
        imirror = imirror + 1;
    end
end
% Plot the connections from stimulated neurons (Stim A)
for iNeuron = DP5.type{1}'
    TargetList = find(DP5.conmat(iNeuron,:) ~= 0);
    TargetList(TargetList == iNeuron) = []; % Do not plot self connections
    mirrori = repmat(DP5.neuronPos(iNeuron,:),[9,1]) + mirroradjust;
    for jNeuron = TargetList
        mirrorj = repmat(DP5.neuronPos(jNeuron,:),[9,1]) + mirroradjust;
        distance = sqrt(sum((repmat(DP5.neuronPos(iNeuron,:),[9,1]) - mirrorj).^2,2));
        minD = find(distance == min(distance),1,'first');
        if minD ~= 5
            % We must connect to a mirror
            
            % Connect from the iNeuron in the main to the jNeuron in the
            % mirror
            line([DP5.neuronPos(iNeuron,1),mirrorj(minD,1)],[DP5.neuronPos(iNeuron,2),mirrorj(minD,2)],'Color',PColor{4},'LineWidth',ConLW)
            
            % Connect from the jNeuron in the main to the iNeuron in the
            % mirror
            minD = 10 - minD;
            line([mirrori(minD,1),DP5.neuronPos(jNeuron,1)],[mirrori(minD,2),DP5.neuronPos(jNeuron,2)],'Color',PColor{4},'LineWidth',ConLW)
        else
            % We can connect without the mirror, so we only need one line
            line([DP5.neuronPos(iNeuron,1),DP5.neuronPos(jNeuron,1)],[DP5.neuronPos(iNeuron,2),DP5.neuronPos(jNeuron,2)],'Color',PColor{4},'LineWidth',ConLW)
        end
    end
end
% Plot the connections from stimulated neurons (Stim B)
for iNeuron = DP5.type{2}'
    TargetList = find(DP5.conmat(iNeuron,:) ~= 0);
    TargetList(TargetList == iNeuron) = []; % Do not plot self connections
    mirrori = repmat(DP5.neuronPos(iNeuron,:),[9,1]) + mirroradjust;
    for jNeuron = TargetList
        mirrorj = repmat(DP5.neuronPos(jNeuron,:),[9,1]) + mirroradjust;
        distance = sqrt(sum((repmat(DP5.neuronPos(iNeuron,:),[9,1]) - mirrorj).^2,2));
        minD = find(distance == min(distance),1,'first');
        if minD ~= 5
            % We must connect to a mirror
            
            % Connect from the iNeuron in the main to the jNeuron in the
            % mirror
            line([DP5.neuronPos(iNeuron,1),mirrorj(minD,1)],[DP5.neuronPos(iNeuron,2),mirrorj(minD,2)],'Color',PColor{4},'LineWidth',ConLW)
            
            % Connect from the jNeuron in the main to the iNeuron in the
            % mirror
            minD = 10 - minD;
            line([mirrori(minD,1),DP5.neuronPos(jNeuron,1)],[mirrori(minD,2),DP5.neuronPos(jNeuron,2)],'Color',PColor{4},'LineWidth',ConLW)
        else
            % We can connect without the mirror, so we only need one line
            line([DP5.neuronPos(iNeuron,1),DP5.neuronPos(jNeuron,1)],[DP5.neuronPos(iNeuron,2),DP5.neuronPos(jNeuron,2)],'Color',PColor{4},'LineWidth',ConLW)
        end
    end
end


% Plot the neuron locations
scatter(DP5.pos{3}(:,1),DP5.pos{3}(:,2),neurDotSz,PColor{1})
scatter(DP5.pos{4}(:,1),DP5.pos{4}(:,2),neurDotSz,PColor{3})
scatter(DP5.pos{1}(:,1),DP5.pos{1}(:,2),neurDotSz,PColor{2},'filled')
scatter(DP5.pos{2}(:,1),DP5.pos{2}(:,2),neurDotSz,PColor{5},'filled')

% Make the legend
xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;
x1 = xMin + 0.4*(xMax - xMin);
x2 = xMin + 0.43*(xMax - xMin);
x3 = xMin + 0.46*(xMax - xMin);
x4 = xMin + 0.64*(xMax - xMin);
y1 = yMin + 0.03*(yMax - yMin);
y2 = yMin + 0.06*(yMax - yMin);
y3 = yMin + 0.105*(yMax - yMin);
y4 = yMin + 0.15*(yMax - yMin);
y5 = yMin + 0.195*(yMax - yMin);
y6 = yMin + 0.225*(yMax - yMin);
fill([x1,x4,x4,x1],[y1,y1,y6,y6],[1,1,1],'EdgeColor','k')
scatter(x2,y5,neurDotSz,PColor{2},'filled')
text(x3,y5,'Stim A','FontSize',LegFS)
scatter(x2,y4,neurDotSz,PColor{5},'filled')
text(x3,y4,'Stim B','FontSize',LegFS)
scatter(x2,y3,neurDotSz,PColor{1},'filled')
text(x3,y3,'Ex','FontSize',LegFS)
scatter(x2,y2,neurDotSz,PColor{3},'filled')
text(x3,y2,'Inh','FontSize',LegFS)

xlabel('X Position','FontSize',AxLabelFS)
ylabel('Y Position','FontSize',AxLabelFS)
set(gca,'FontSize',UnitFS)
set(gca,'Layer','Top')
title('Example Network','FontSize',TitleFS)

xlim([xMin,xMax])
ylim([yMin,yMax])



% Make the example spike rasters

for iStim = 1:3
    
    if iStim == 1
        iCol = 2;
    elseif iStim == 2
        iCol = 1;
    elseif iStim == 3
        iCol = 3;
    end
    
    f = subplot('Position',[LeftCoord1(1 + iCol),BottomCoord1(1),width3,height1]);
    hold on
    
    % Set limits
    xMin = -20;
    xMax = 60;
    yMin = 0;
    yMax = 1000;
    
    
    line([0,0],[0,1000],'Color','k','LineStyle','--','LineWidth',StimTLW)
    
    % Plot the spikes
    scatter(DP6.spikes{1,iStim}(:,2),DP6.spikes{1,iStim}(:,1),spikeDotSz,'filled','MarkerFaceColor',PColor{2})
    scatter(DP6.spikes{2,iStim}(:,2),DP6.spikes{2,iStim}(:,1),spikeDotSz,'filled','MarkerFaceColor',PColor{5})
    scatter(DP6.spikes{3,iStim}(:,2),DP6.spikes{3,iStim}(:,1),spikeDotSz,'filled','MarkerFaceColor',PColor{1})
    scatter(DP6.spikes{4,iStim}(:,2),DP6.spikes{4,iStim}(:,1),spikeDotSz,'filled','MarkerFaceColor',PColor{3})
    
    % Make the legend
    if iCol == 1
        x1 = xMin + 0.65*(xMax - xMin);
        x2 = xMin + 0.69*(xMax - xMin);
        x3 = xMin + 0.73*(xMax - xMin);
        x4 = xMin + 0.95*(xMax - xMin);
        y1 = yMin + 0.03*(yMax - yMin);
        y2 = yMin + 0.06*(yMax - yMin);
        y3 = yMin + 0.105*(yMax - yMin);
        y4 = yMin + 0.15*(yMax - yMin);
        y5 = yMin + 0.195*(yMax - yMin);
        y6 = yMin + 0.225*(yMax - yMin);
        fill([x1,x4,x4,x1],[y1,y1,y6,y6],[1,1,1],'EdgeColor','k')
        scatter(x2,y5,neurDotSz,PColor{2},'filled')
        text(x3,y5,'Stim A','FontSize',LegFS)
        scatter(x2,y4,neurDotSz,PColor{5},'filled')
        text(x3,y4,'Stim B','FontSize',LegFS)
        scatter(x2,y3,neurDotSz,PColor{1},'filled')
        text(x3,y3,'Ex','FontSize',LegFS)
        scatter(x2,y2,neurDotSz,PColor{3},'filled')
        text(x3,y2,'Inh','FontSize',LegFS)
    end
    
    xlim([xMin,xMax])
    ylim([yMin,yMax])
    set(gca,'FontSize',UnitFS)
    
    xlabel('Time Rel. to Stim (ms)','FontSize',AxLabelFS)
    if iCol == 1
        title('Stimulus A','FontSize',TitleFS)
        ylabel('Neuron Number','FontSize',AxLabelFS)
    elseif iCol == 2
        title(['Example Spike Raster',char(10),'Stimulus B'],'FontSize',TitleFS)
        set(gca,'YTickLabel',[])
    elseif iCol == 3
        title('Stimulus A and B','FontSize',TitleFS)
        set(gca,'YTickLabel',[])
    end
    
end


% Make the PID heat maps

% Set the minimum value reset for plotting in log scale
minReset = 10^-5;
minNonZero = min(DP7.PIDHeatMap(DP7.PIDHeatMap(:) > (100*eps))); % Look for the smallest non-zero result, not counting rounding errors
minReset = max([minReset,minNonZero]);
DP7.PIDHeatMap(DP7.PIDHeatMap < minReset) = minReset;


for iPID = 1:4
    
    % Find the color limits
    cMax = -inf;
    cMin = inf;
    for iPlot = 1:nHM
        cMax = max([cMax,max(max(max(squeeze(DP7.PIDHeatMap(:,:,HMTBins(iPlot),iPID)))))]);
        cMin = min([cMin,min(min(min(squeeze(DP7.PIDHeatMap(:,:,HMTBins(iPlot),iPID)))))]);
    end
    
    
    for iPlot = 1:nHM
        f = subplot('Position',[LeftCoord2(iPlot),BottomCoord2(iPID),width2,height2]);
        hold on
        
        imagesc(log10(squeeze(DP7.PIDHeatMap(:,:,HMTBins(iPlot),iPID))'),log10([cMin,cMax]))
        
        line([(0.25*DP7.nBins) + 0.5,(0.25*DP7.nBins) + 0.5],[0.5,DP7.nBins + 0.5],'Color','k','LineWidth',NeuronLineLW)
        line([(0.25*DP7.nBins) + 0.5,(0.25*DP7.nBins) + 0.5],[0.5,DP7.nBins + 0.5],'Color','w','LineStyle','--','LineWidth',NeuronLineLW)
        line([(0.75*DP7.nBins) + 0.5,(0.75*DP7.nBins) + 0.5],[0.5,DP7.nBins + 0.5],'Color','k','LineWidth',NeuronLineLW)
        line([(0.75*DP7.nBins) + 0.5,(0.75*DP7.nBins) + 0.5],[0.5,DP7.nBins + 0.5],'Color','w','LineStyle','--','LineWidth',NeuronLineLW)
        
        xlim([0.5,DP7.nBins + 0.5])
        ylim([0.5,DP7.nBins + 0.5])
        
        set(gca,'FontSize',UnitFS)
        if iPID == 1
            if iPlot > 1
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            end
            if iPlot == 1
                set(gca,'XTick',[0.5,(DP7.nBins/2) + 0.5,DP7.nBins + 0.5])
                set(gca,'XTickLabel',{'0','0.5','1'})
                set(gca,'YTick',[0.5,(DP7.nBins/2) + 0.5,DP7.nBins + 0.5])
                set(gca,'YTickLabel',{'0','0.5','1'})
                xlabel('X Position','FontSize',AxLabelFS)
                ylabel('Y Position','FontSize',AxLabelFS)
            end
            if iPlot == round(nHM/2)
                title(['Average Stimulus Encoding (',PIDNamesLong{iPID},' All Models)',char(10),...
                    num2str(DP7.Time(HMTBins(iPlot),1)),' to ',num2str(DP7.Time(HMTBins(iPlot),2)),' ms'],'FontSize',TitleFS)
            else
                title([num2str(DP7.Time(HMTBins(iPlot),1)),' to ',num2str(DP7.Time(HMTBins(iPlot),2)),' ms'],'FontSize',TitleFS)
            end
        else
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            if iPlot == round(nHM/2)
                title(['Average Stimulus Encoding (',PIDNamesLong{iPID},' All Models)'],'FontSize',TitleFS)
            end
        end
        
        
        
    end
    
    % Add the colorbar
    yticks = [];
    precision = 1;
    while length(yticks) < 2
        yticks = (ceil(log10(cMin)/precision):floor(log10(cMax)/precision))*precision;
        precision = precision/2;
    end
    yticklabels = cell([1,length(yticks)]);
    for iTick = 1:length(yticks)
        yticklabels{iTick} = num2str(yticks(iTick));
    end
    hbar = colorbar;
    set(hbar,'location','manual',...
        'position',[LeftCoord2(nHM + 1),BottomCoord2(iPID),cbwidth,height2],...
        'fontsize',CBFS,...
        'YTick',yticks,...
        'YTickLabel',yticklabels,...
        'ylim',log10([cMin,cMax]));
    ylabel(hbar,['log_{10} ',PIDNamesShort{iPID},' (bits)'],'FontSize',AxLabelFS)
    
end



% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo3Panel2')
cd(SimDir)
