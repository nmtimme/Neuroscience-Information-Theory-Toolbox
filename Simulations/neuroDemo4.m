%% neuroDemo4
%
%   neuroDemo4 reproduces the small network simulations. Note, the scratch 
%   directory in the simulations folder will be used to save data. 
%   These simulations will take up about 1.5 hours to run and it will 
%   require about 60 GB of space.
%
% Other m-files required: modeltype4, instinfo, inverseFNoise,
% simpleModels_FSI, simpleModel_RS
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% November 2016; Last revision: 1-Nov-2016

%% Set Parameters

disp('Note that this simulation will take about 1.5 hours to run and it requires about 60 GB of space.')

% Set the number of models to run for each subtype
nModelsPerSubtype = 20;

% Set the number of correlation values to scan for the correlated input
% subtype
nScan = 11;

% Set the pulse size
pulseSize = 500;

% Set the parameters for the various subtypes
nSTs = 4 + nScan;
weights = zeros([nSTs,5]);
neuronTypes = cell([nSTs,1]);
backgroundCurrent = zeros([nSTs,1]);
varStim = zeros([nSTs,4]);

%%%%%%%%%%%%%%%%%%%%%%%
% StimA Gate

% Set the weights
weights(1,:) = [200,0,0,0,0];

% Set the neuron types
neuronTypes{1} = {'Ex','Ex','Ex','Ex'};

% Set the background current
backgroundCurrent(1) = 0;

% Set the variable stimulation
varStim(1,:) = 0.25*ones([1,4]);

%%%%%%%%%%%%%%%%%%%%%%%
% NOR Gate

% Set the weights
weights(2,:) = [-30,-30,0,0,0];

% Set the neuron types
neuronTypes{2} = {'Inh','Inh','Ex','Ex'};

% Set the background current
backgroundCurrent(2) = 200;

% Set the variable stimulation
varStim(2,:) = 0.25*ones([1,4]);

%%%%%%%%%%%%%%%%%%%%%%%
% NOR Gate w/ No Background

% Set the weights
weights(3,:) = [-30,-30,0,0,0];

% Set the neuron types
neuronTypes{3} = {'Inh','Inh','Ex','Ex'};

% Set the background current
backgroundCurrent(3) = 0;

% Set the variable stimulation
varStim(3,:) = 0.25*ones([1,4]);

%%%%%%%%%%%%%%%%%%%%%%%
% XOR Gate

% Set the weights
weights(4,:) = [200,200,20,20,-500]; % Note, changed 50 to 20, changed -250 to -500

% Set the neuron types
neuronTypes{4} = {'Ex','Ex','Inh','Ex'};

% Set the background current
backgroundCurrent(4) = -100;

% Set the variable stimulation
varStim(4,:) = 0.25*ones([1,4]);

%%%%%%%%%%%%%%%%%%%%%%%
% OR Gate w/ Variable Correlation

up = linspace(0,0.5,nScan);
down = linspace(0.5,0,nScan);
for iCor = 1:nScan
    % Set the weights
    weights(4 + iCor,:) = [100,100,0,0,0];
    
    % Set the neuron types
    neuronTypes{4 + iCor} = {'Ex','Ex','Ex','Ex'};
    
    % Set the background current
    backgroundCurrent(4 + iCor) = 0;
    
    % Set the variable stimulation
    varStim(4 + iCor,:) = [up(iCor),down(iCor),down(iCor),up(iCor)];
end

% Set the bin size in milliseconds
BinSize = 25; 

% Set the time window around the start of the stimuli to analyze in
% milliseconds
MaxLead = 200;
MaxLag = 700;

% Get the scratch directory to use to store data
if ispc
    ScratchDir = [pwd,'\Scratch\ModelType4\'];
elseif isunix
    ScratchDir = [pwd,'/Scratch/ModelType4/'];
elseif ismac
    ScratchDir = [pwd,'/Scratch/ModelType4/'];
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
            
        [mV,spkTimes,memNoise,current,stim,modelParams] = modeltype4(pulseSize, weights(iSubType,:), neuronTypes{iSubType}, backgroundCurrent(iSubType), ...
            'varStim', varStim(iSubType,:));
        
        save([ScratchDir,'ModelType4ST',num2str(iSubType),'Run',num2str(iModel),'.mat'],'mV','spkTimes','memNoise','current','stim','modelParams')
        
        disp(['Finished Generating Data for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Convert the Raw Data to DataRaster Format

% Go to the directory with the analysis software
cd(AnaDir)

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType4ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])
        
        % Find the number of neurons
        nNeurons = size(spkTimes,1);
        
        % Make the Data Raster
        DataRaster = cell([2,1]);
        
        % Error check that the stimulation times are the same
        if ~isequal(stim.times{1},stim.times{2})
            error('Stimulation times are not identical.')
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
            DataRaster{2}(i,:,:) = reshape(stim.types{i},[1,1,length(stim.types{i})]);
        end
        
        % Record the binsize
        modelParams.binsize = BinSize;
        modelParams.maxlead = MaxLead;
        modelParams.maxlag = MaxLag;
        modelParams.Time = timeboundaries;
        
        % Save the Data
        save([ScratchDir,'ModelType4ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'],'DataRaster','modelParams')
        
        disp(['Finished Converting Data to DataRaster Format for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Convert the DataRaster Format Data to StatesRasters

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType4ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])
        
        % Make the stating method assignment array (making sure to count the
        % number of pulses that were used)
        nPulses = length(modelParams.pulseSize);
        pSet = (1/nPulses)*ones([1,nPulses]);
        MethodAssign = {2,1,'Nat',{};...
            2,2,'Nat',{};...
            1,1,'UniCB',{nPulses};...
            1,2,'UniCB',{nPulses};...
            1,3,'UniCB',{nPulses};...
            1,4,'UniCB',{nPulses}};
        
        % Convert the data
        StatesRaster = data2states(DataRaster,MethodAssign);
        
        % Save the Data
        save([ScratchDir,'ModelType4ST',num2str(iSubType),'Run',num2str(iModel),'SR.mat'],'StatesRaster','modelParams')
        
        disp(['Finished Converting Data to StatesRaster Format for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Perform the Information Theory Analysis

InfoResults = cell([nSTs,nModelsPerSubtype]);
SubTypeParams = cell([nSTs,1]);
for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType4ST',num2str(iSubType),'Run',num2str(iModel),'SR.mat'])
        
        % Save the subtype parameters. We assume all runs of a subtype used the
        % same parameters.
        if iModel == 1
            SubTypeParams{iSubType} = modelParams;
        end
        
        % Figure out how many time steps we have in this data set
        nT = size(StatesRaster{1},2);
        
        % Preallocate variables
        PIDVals = NaN([4,nT]);
        
        % Stimulus encoding by output neuron
        Method = '2PID';
        for iT = 1:nT
            VariableIDs = {1,4,iT;2,1,1;2,2,1};
            PIDVals(:,iT) = instinfo(StatesRaster, Method, VariableIDs);
        end
        
        InfoResults{iSubType,iModel} = PIDVals;
        
        disp(['Finished Analyzing for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

% Save the Results
save([ScratchDir,'ModelType4InfoResults.mat'],'InfoResults','SubTypeParams')

%% Gather Data for Plotting

% Get Example Voltage Traces

% Select the subtypes of interest
targetSTs = [1,2,3,4,10];

DP1 = struct;

DP1.spkTimes = cell([length(targetSTs),4]);
DP1.time = cell([length(targetSTs),1]);
DP1.mV = cell([length(targetSTs),4]);
DP1.current = cell([length(targetSTs),2]);
for iSubType = 1:length(targetSTs)
    
    % Load the first run
    load([ScratchDir,'ModelType4ST',num2str(targetSTs(iSubType)),'Run1.mat'])
    
    % Get the total time vector
    time = (1:size(mV,2))*modelParams.tau;
    
    % Set the minimum and maximum times based on the subtype
    stimset = 10;
    
    ExampleTime = length(modelParams.pulseSize)*modelParams.intrastimT + (length(modelParams.pulseSize) + 1)*modelParams.interstimT;
    OneSetTime = length(modelParams.pulseSize)*modelParams.intrastimT + length(modelParams.pulseSize)*modelParams.interstimT;
    MinTime = (stimset - 1)*OneSetTime;
    MaxTime = MinTime + ExampleTime;
    
    % Get the time vector
    DP1.time{iSubType} = time((time > MinTime) & (time <= MaxTime));
    
    for i = 1:4
        
        % Get the spike times,
        DP1.spkTimes{iSubType,i} = spkTimes{i}(ismember(spkTimes{i},DP1.time{iSubType}));
        
        % Get the membrane voltage
        DP1.mV{iSubType,i} = mV(i,(time > MinTime) & (time <= MaxTime));
        
    end
    
    % Get the stimulus current
    for iStim = 1:2
        DP1.current{iSubType,iStim} = stim.currentOriginal(iStim,(time > MinTime) & (time <= MaxTime));
    end
    
    % Correct offset in time for plotting
    DP1.time{iSubType} = DP1.time{iSubType} - MinTime;
    for i = 1:4
        DP1.spkTimes{iSubType,i} = DP1.spkTimes{iSubType,i} - MinTime;
    end
end



% Get the information results through time for the StimA gate

DP2 = struct;

iSubType = 1;

modelParams = SubTypeParams{iSubType};

% Get the time values
DP2.time = mean(modelParams.Time,2)';
nBins = length(DP2.time);

% Get all the PID information values and p-values
AllInfo = zeros([nModelsPerSubtype,nBins,4]);
for iModel = 1:nModelsPerSubtype
    for iPID = 1:4
        AllInfo(iModel,:,iPID) = InfoResults{iSubType,iModel}(iPID,:);
    end
end

DP2.boxpoints = NaN([5,nBins,3]);
DP2.meanerror = NaN([5,nBins,3]);
for iBin = 1:nBins
    for iPID = 1:4
        temp = sort(AllInfo(:,iBin,iPID),'descend');
        DP2.boxpoints(1,iBin,iPID) = temp(1);
        DP2.boxpoints(2,iBin,iPID) = temp(floor(nModelsPerSubtype*0.75));
        DP2.boxpoints(3,iBin,iPID) = median(temp);
        DP2.boxpoints(4,iBin,iPID) = temp(ceil(nModelsPerSubtype*0.25));
        DP2.boxpoints(5,iBin,iPID) = temp(end);
        DP2.meanerror(1,iBin,iPID) = mean(temp);
        DP2.meanerror(2,iBin,iPID) = mean(temp) + std(temp)/sqrt(length(temp));
        DP2.meanerror(3,iBin,iPID) = mean(temp) - std(temp)/sqrt(length(temp));
        DP2.meanerror(4,iBin,iPID) = mean(temp) + std(temp);
        DP2.meanerror(5,iBin,iPID) = mean(temp) - std(temp);
    end
end


% Get the information results through time for the NOR gate

DP3 = struct;
DP3.time = cell([2,1]);
DP3.boxpoints = cell([2,1]);
DP3.meanerror = cell([2,1]);
DP3.nsig = cell([2,1]);

for iSubType = 2:3
    
    modelParams = SubTypeParams{iSubType};
    
    % Get the time values
    DP3.time{iSubType - 1} = mean(modelParams.Time,2)';
    nBins = length(DP3.time{iSubType - 1});
    
    % Get all the PID information values and p-values
    AllInfo = zeros([nModelsPerSubtype,nBins,4]);
    for iModel = 1:nModelsPerSubtype
        for iPID = 1:4
            AllInfo(iModel,:,iPID) = InfoResults{iSubType,iModel}(iPID,:);
        end
    end
    
    DP3.boxpoints{iSubType - 1} = NaN([5,nBins,3]);
    DP3.meanerror{iSubType - 1} = NaN([5,nBins,3]);
    for iBin = 1:nBins
        for iPID = 1:4
            temp = sort(AllInfo(:,iBin,iPID),'descend');
            DP3.boxpoints{iSubType - 1}(1,iBin,iPID) = temp(1);
            DP3.boxpoints{iSubType - 1}(2,iBin,iPID) = temp(floor(nModelsPerSubtype*0.75));
            DP3.boxpoints{iSubType - 1}(3,iBin,iPID) = median(temp);
            DP3.boxpoints{iSubType - 1}(4,iBin,iPID) = temp(ceil(nModelsPerSubtype*0.25));
            DP3.boxpoints{iSubType - 1}(5,iBin,iPID) = temp(end);
            DP3.meanerror{iSubType - 1}(1,iBin,iPID) = mean(temp);
            DP3.meanerror{iSubType - 1}(2,iBin,iPID) = mean(temp) + std(temp)/sqrt(length(temp));
            DP3.meanerror{iSubType - 1}(3,iBin,iPID) = mean(temp) - std(temp)/sqrt(length(temp));
            DP3.meanerror{iSubType - 1}(4,iBin,iPID) = mean(temp) + std(temp);
            DP3.meanerror{iSubType - 1}(5,iBin,iPID) = mean(temp) - std(temp);
        end
    end
end


% Get the information results through time for the XOR gate

DP4 = struct;

iSubType = 4;

modelParams = SubTypeParams{iSubType};

% Get the time values
DP4.time = mean(modelParams.Time,2)';
nBins = length(DP4.time);

% Get all the PID information values and p-values
AllInfo = zeros([nModelsPerSubtype,nBins,4]);
for iModel = 1:nModelsPerSubtype
    for iPID = 1:4
        AllInfo(iModel,:,iPID) = InfoResults{iSubType,iModel}(iPID,:);
    end
end

DP4.boxpoints = NaN([5,nBins,3]);
DP4.meanerror = NaN([5,nBins,3]);
for iBin = 1:nBins
    for iPID = 1:4
        temp = sort(AllInfo(:,iBin,iPID),'descend');
        DP4.boxpoints(1,iBin,iPID) = temp(1);
        DP4.boxpoints(2,iBin,iPID) = temp(floor(nModelsPerSubtype*0.75));
        DP4.boxpoints(3,iBin,iPID) = median(temp);
        DP4.boxpoints(4,iBin,iPID) = temp(ceil(nModelsPerSubtype*0.25));
        DP4.boxpoints(5,iBin,iPID) = temp(end);
        DP4.meanerror(1,iBin,iPID) = mean(temp);
        DP4.meanerror(2,iBin,iPID) = mean(temp) + std(temp)/sqrt(length(temp));
        DP4.meanerror(3,iBin,iPID) = mean(temp) - std(temp)/sqrt(length(temp));
        DP4.meanerror(4,iBin,iPID) = mean(temp) + std(temp);
        DP4.meanerror(5,iBin,iPID) = mean(temp) - std(temp);
    end
end


% Get the information results through time for the OR gate for specific correlation values

% Set the desired subtypes
targetSTs = [5,10,15];
ntargetSTs = length(targetSTs);

DP5 = struct;
DP5.time = cell([ntargetSTs,1]);
DP5.boxpoints = cell([ntargetSTs,1]);
DP5.meanerror = cell([ntargetSTs,1]);
DP5.nsig = cell([ntargetSTs,1]);

for iSubType = 1:ntargetSTs
    
    modelParams = SubTypeParams{targetSTs(iSubType)};
    
    % Get the time values
    DP5.time{iSubType} = mean(modelParams.Time,2)';
    nBins = length(DP5.time{iSubType});
    
    % Get all the PID information values and p-values
    AllInfo = zeros([nModelsPerSubtype,nBins,4]);
    for iModel = 1:nModelsPerSubtype
        for iPID = 1:4
            AllInfo(iModel,:,iPID) = InfoResults{targetSTs(iSubType),iModel}(iPID,:);
        end
    end
    
    DP5.boxpoints{iSubType} = NaN([5,nBins,3]);
    DP5.meanerror{iSubType} = NaN([5,nBins,3]);
    for iBin = 1:nBins
        for iPID = 1:4
            temp = sort(AllInfo(:,iBin,iPID),'descend');
            DP5.boxpoints{iSubType}(1,iBin,iPID) = temp(1);
            DP5.boxpoints{iSubType}(2,iBin,iPID) = temp(floor(nModelsPerSubtype*0.75));
            DP5.boxpoints{iSubType}(3,iBin,iPID) = median(temp);
            DP5.boxpoints{iSubType}(4,iBin,iPID) = temp(ceil(nModelsPerSubtype*0.25));
            DP5.boxpoints{iSubType}(5,iBin,iPID) = temp(end);
            DP5.meanerror{iSubType}(1,iBin,iPID) = mean(temp);
            DP5.meanerror{iSubType}(2,iBin,iPID) = mean(temp) + std(temp)/sqrt(length(temp));
            DP5.meanerror{iSubType}(3,iBin,iPID) = mean(temp) - std(temp)/sqrt(length(temp));
            DP5.meanerror{iSubType}(4,iBin,iPID) = mean(temp) + std(temp);
            DP5.meanerror{iSubType}(5,iBin,iPID) = mean(temp) - std(temp);
        end
    end
end



% Get the information results for the OR gate across correlation values

DP6 = struct;

DP6.time = cell([nSTs - 4,1]);
DP6.boxpoints = cell([nSTs - 4,1]);
DP6.nsig = cell([nSTs - 4,1]);
DP6.meanerror = cell([nSTs - 4,1]);
DP6.InhibVar = NaN([10,nSTs - 4,3]);
DP6.CorParam = NaN([1,nSTs - 4]);
for iSubType = 5:nSTs
    
    modelParams = SubTypeParams{iSubType};
    
    
    % Get the time values
    DP6.time{iSubType - 4} = mean(modelParams.Time,2)';
    nBins = length(DP6.time{iSubType - 4});
    
    % Get all the PID information values and p-values
    AllInfo = zeros([nModelsPerSubtype,nBins,4]);
    for iModel = 1:nModelsPerSubtype
        for iPID = 1:4
            AllInfo(iModel,:,iPID) = InfoResults{iSubType,iModel}(iPID,:);
        end
    end
    
    % Get average information during stimulus as a function of inhibition
    for iPID = 1:4
        temp = AllInfo(:,(DP6.time{iSubType - 4} > 0) & (DP6.time{iSubType - 4} < modelParams.intrastimT),iPID);
        temp = sort(temp(:),'descend');
        DP6.InhibVar(1,iSubType - 4,iPID) = temp(1);
        DP6.InhibVar(2,iSubType - 4,iPID) = temp(floor(nModelsPerSubtype*0.75));
        DP6.InhibVar(3,iSubType - 4,iPID) = median(temp);
        DP6.InhibVar(4,iSubType - 4,iPID) = temp(ceil(nModelsPerSubtype*0.25));
        DP6.InhibVar(5,iSubType - 4,iPID) = temp(end);
        DP6.InhibVar(6,iSubType - 4,iPID) = mean(temp);
        DP6.InhibVar(7,iSubType - 4,iPID) = mean(temp) + std(temp)/sqrt(length(temp));
        DP6.InhibVar(8,iSubType - 4,iPID) = mean(temp) - std(temp)/sqrt(length(temp));
        DP6.InhibVar(9,iSubType - 4,iPID) = mean(temp) + std(temp);
        DP6.InhibVar(10,iSubType - 4,iPID) = mean(temp) - std(temp);
    end
    
    % Record the inhibitory strength
    DP6.CorParam(iSubType - 4) = modelParams.varStim(1) - 0.25;
    
    
    
end


% Save the results
save([ScratchDir,'ModelType4PlotData'],'DP1','DP2','DP3','DP4','DP5','DP6')



%% Make the figure for the StimA gate

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.18;

% Set the top and bottom margins in inches
tmargin = 0.23;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace = [0.4,0.4,0.4];

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set Paper size
papersize = [7.4 3.5];

% Set the figure panel dimensions (number of horizontal figures by number
% of vertical figures)
figdim = [4,2];

% Set the height proportions for each row
heightprop = [1,1];

% Set the width proportions for each column
widthprop = [1,1,1,1];

% calculate the width and height based on dimensions above in inches
widthtot = papersize(1) - lmargin - rmargin - sum(hspace);
widthunit = widthtot/sum(widthprop);
width = widthunit*widthprop;
heighttot = papersize(2) - tmargin - bmargin - (figdim(2) - 1)*vspace;
heightunit = heighttot/sum(heightprop);
height = heightunit*heightprop;

LeftCoord = cumsum([0,width(1:(end - 1)) + hspace]) + lmargin;
BottomCoord = cumsum([0,fliplr(height + vspace)]) + bmargin;
BottomCoord(end) = [];
BottomCoord = fliplr(BottomCoord);

% Convert to fraction of page sizes
LeftCoord = LeftCoord/papersize(1);
BottomCoord = BottomCoord/papersize(2);
width = width/papersize(1);
height = height/papersize(2);
hspace = hspace/papersize(1);

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
StimLabelFS = 6;

% Set the linewidths
VoltageLW = 0.6;
BPThickLW = 3;
BPThinLW = 0.6;
StimLineLW = 1;
DelayLW = 2;
AxLW = 0.4;
SpikeLW = 0.4;

% Set dot sizes
BPDotSz = 8;

% Make PID names
PIDNamesLong = {'Redundancy','Unique A','Unique B','Synergy'};
PIDNamesShort = {'Red','Uniq A','Uniq B','Syn'};

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');

% Make the spike rasters

iST = 1;

f = subplot('Position',[LeftCoord(3),BottomCoord(1),width(3) + width(4) + hspace(3),height(1)]);
hold on

% Set the limits
xMin = 0;
xMax = max(DP1.time{iST});
yMax = max([DP1.current{iST,1},DP1.current{iST,2}]);
yMin = min([DP1.current{iST,1},DP1.current{iST,2}]);
dy = yMax - yMin;

% Plot the stimulus current
plot(DP1.time{iST},DP1.current{iST,1},'k','LineWidth',VoltageLW)
plot(DP1.time{iST},DP1.current{iST,2} + yMax + dy,'k','LineWidth',VoltageLW) % Assume inputs have same magnitude

text(xMin + 0.18*(xMax - xMin),0.5*dy,'Stimulus A','FontSize',StimLabelFS)
text(xMin + 0.18*(xMax - xMin),2.5*dy,'Stimulus B','FontSize',StimLabelFS)

line([xMin + 0.1*(xMax - xMin),xMin + 0.1*(xMax - xMin)],[0.5*dy,1.5*dy],'Color','k','LineWidth',VoltageLW)
text(xMin + 0.09*(xMax - xMin),dy,[num2str(dy),' pA'],'FontSize',LegFS,'HorizontalAlignment','right')

yMax = yMax + dy + yMax;

% Plot the spikes
for iSpike = 1:length(DP1.spkTimes{iST,1})
    line([DP1.spkTimes{iST,1}(iSpike),DP1.spkTimes{iST,1}(iSpike)],[yMax + dy,yMax + 2*dy],'Color',PColor{2},'LineWidth',SpikeLW)
end
for iSpike = 1:length(DP1.spkTimes{iST,2})
    line([DP1.spkTimes{iST,2}(iSpike),DP1.spkTimes{iST,2}(iSpike)],[yMax + 3*dy,yMax + 4*dy],'Color',PColor{3},'LineWidth',SpikeLW)
end
for iSpike = 1:length(DP1.spkTimes{iST,4})
    line([DP1.spkTimes{iST,4}(iSpike),DP1.spkTimes{iST,4}(iSpike)],[yMax + 5*dy,yMax + 6*dy],'Color',PColor{1},'LineWidth',SpikeLW)
end

text(xMin + 0.18*(xMax - xMin),yMax + 2.5*dy,'Neuron E1 Spikes','FontSize',StimLabelFS,'Color',PColor{2})
text(xMin + 0.18*(xMax - xMin),yMax + 4.5*dy,'Neuron E2 Spikes','FontSize',StimLabelFS,'Color',PColor{3})
text(xMin + 0.18*(xMax - xMin),yMax + 6.5*dy,'Neuron E3 Spikes','FontSize',StimLabelFS,'Color',PColor{1})


xlim([xMin,xMax])
ylim([yMin - dy,yMax + 7*dy])
set(gca,'FontSize',UnitFS)
xlabel('Time (ms)','FontSize',AxLabelFS)
title('Example Stimuli and Spike Rasters','FontSize',TitleFS)
set(gca,'ytick',[])


% Plot the PID through Time Results

yMax = -inf;
yMin = inf;
for iPID = 1:4
    yMax = max([max(max(DP2.meanerror(:,:,iPID))),yMax]);
    yMin = min([min(min(DP2.meanerror(:,:,iPID))),yMin]);
end

dy = yMax - yMin;
yMax = yMax + 0.1*dy;
yMin = yMin - 0.1*dy;

xMin = DP2.time(1) - 0.5*(DP2.time(2) - DP2.time(1));
xMax = DP2.time(end) + 0.5*(DP2.time(2) - DP2.time(1));

for iPID = 1:4

    f = subplot('Position',[LeftCoord(iPID),BottomCoord(2),width(iPID),height(2)]);
    hold on
    
    nBins = size(DP2.meanerror,2);
    
    line([0,0],[yMin,yMax],'Color',PColor{3},'LineStyle','--','LineWidth',StimLineLW)
    text(-15,mean([yMin,yMax]),'Stim On','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',90,'FontSize',StimLabelFS,'Color',PColor{3})
    line([500,500],[yMin,yMax],'Color',PColor{2},'LineStyle','--','LineWidth',StimLineLW)
    text(500 + 15,mean([yMin,yMax]),'Stim Off','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-90,'FontSize',StimLabelFS,'Color',PColor{2})
    
    for iBin = 1:nBins
        line([DP2.time(iBin),DP2.time(iBin)],[DP2.meanerror(4,iBin,iPID),DP2.meanerror(5,iBin,iPID)],'Color','k','LineWidth',BPThinLW)
        scatter(DP2.time(iBin),DP2.meanerror(1,iBin,iPID),BPDotSz,'k','filled')
    end
    
    ylim([yMin,yMax])
    xlim([xMin,xMax])
    
    
    set(gca,'FontSize',UnitFS)
    xlabel('Time (ms)','FontSize',AxLabelFS)
    ylabel([PIDNamesShort{iPID},' (bits)'],'FontSize',AxLabelFS)
    title(PIDNamesLong{iPID},'FontSize',TitleFS)
    
end


% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo4StimA')
cd(SimDir)



%% Make the figure for the NOR gate

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.18;

% Set the top and bottom margins in inches
tmargin = 0.23;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace = [0.4,0.4,0.4];

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set Paper size
papersize = [7.4 4.75];

% Set the figure panel dimensions (number of horizontal figures by number
% of vertical figures)
figdim = [4,3];

% Set the height proportions for each row
heightprop = [1,1,1];

% Set the width proportions for each column
widthprop = [1,1,1,1];

% calculate the width and height based on dimensions above in inches
widthtot = papersize(1) - lmargin - rmargin - sum(hspace);
widthunit = widthtot/sum(widthprop);
width = widthunit*widthprop;
heighttot = papersize(2) - tmargin - bmargin - (figdim(2) - 1)*vspace;
heightunit = heighttot/sum(heightprop);
height = heightunit*heightprop;

LeftCoord = cumsum([0,width(1:(end - 1)) + hspace]) + lmargin;
BottomCoord = cumsum([0,fliplr(height + vspace)]) + bmargin;
BottomCoord(end) = [];
BottomCoord = fliplr(BottomCoord);

% Convert to fraction of page sizes
LeftCoord = LeftCoord/papersize(1);
BottomCoord = BottomCoord/papersize(2);
width = width/papersize(1);
height = height/papersize(2);
hspace = hspace/papersize(1);

% Set Colors
PColor = {[0,0,1],[1,0,0],[0,210/255,50/255],[0,0,1],[150/255,0,200/255],[200/255,0,255/255],[0/255,255/255,185/255],[238/255,130/255,238/255]};
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

% Set the linewidths
VoltageLW = 0.6;
BPThickLW = 3;
BPThinLW = 0.6;
StimLineLW = 1;
DelayLW = 2;
AxLW = 0.4;
SpikeLW = 0.1;

% Set dot sizes
BPDotSz = 8;

% Make PID names
PIDNamesLong = {'Redundancy','Unique A','Unique B','Synergy'};
PIDNamesShort = {'Red','Uniq A','Uniq B','Syn'};

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');

% Make the spike rasters

for iST = 2:3
    
    f = subplot('Position',[LeftCoord(3),BottomCoord(iST - 1),width(3) + width(4) + hspace(3),height(iST - 1)]);
    hold on
    
    % Set the limits
    xMin = 0;
    xMax = max(DP1.time{iST});
    yMax = max([DP1.current{iST,1},DP1.current{iST,2}]);
    yMin = min([DP1.current{iST,1},DP1.current{iST,2}]);
    dy = yMax - yMin;
    
    % Plot the stimulus current
    plot(DP1.time{iST},DP1.current{iST,1},'k','LineWidth',VoltageLW)
    plot(DP1.time{iST},DP1.current{iST,2} + yMax + dy,'k','LineWidth',VoltageLW) % Assume inputs have same magnitude
    
    text(xMin + 0.18*(xMax - xMin),0.5*dy,'Stimulus A','FontSize',StimLabelFS)
    text(xMin + 0.18*(xMax - xMin),2.5*dy,'Stimulus B','FontSize',StimLabelFS)
    
    line([xMin + 0.1*(xMax - xMin),xMin + 0.1*(xMax - xMin)],[0.5*dy,1.5*dy],'Color','k','LineWidth',VoltageLW)
    text(xMin + 0.09*(xMax - xMin),dy,[num2str(dy),' pA'],'FontSize',LegFS,'HorizontalAlignment','right')
    
    yMax = yMax + dy + yMax;
    
    % Plot the spikes
    for iSpike = 1:length(DP1.spkTimes{iST,1})
        line([DP1.spkTimes{iST,1}(iSpike),DP1.spkTimes{iST,1}(iSpike)],[yMax + dy,yMax + 2*dy],'Color',PColor{2},'LineWidth',SpikeLW)
    end
    for iSpike = 1:length(DP1.spkTimes{iST,2})
        line([DP1.spkTimes{iST,2}(iSpike),DP1.spkTimes{iST,2}(iSpike)],[yMax + 3*dy,yMax + 4*dy],'Color',PColor{3},'LineWidth',SpikeLW)
    end
    for iSpike = 1:length(DP1.spkTimes{iST,4})
        line([DP1.spkTimes{iST,4}(iSpike),DP1.spkTimes{iST,4}(iSpike)],[yMax + 5*dy,yMax + 6*dy],'Color',PColor{1},'LineWidth',SpikeLW)
    end
    
    text(xMin + 0.18*(xMax - xMin),yMax + 2.5*dy,'Neuron I1 Spikes','FontSize',StimLabelFS,'Color',PColor{2})
    text(xMin + 0.18*(xMax - xMin),yMax + 4.5*dy,'Neuron I2 Spikes','FontSize',StimLabelFS,'Color',PColor{3})
    text(xMin + 0.18*(xMax - xMin),yMax + 6.5*dy,'Neuron E1 Spikes','FontSize',StimLabelFS,'Color',PColor{1})
    
    
    xlim([xMin,xMax])
    ylim([yMin - dy,yMax + 7*dy])
    set(gca,'FontSize',UnitFS)
    xlabel('Time (ms)','FontSize',AxLabelFS)
    set(gca,'ytick',[])
    
    if iST == 2
        title('Example Stimuli and Spike Rasters (Background ON)','FontSize',TitleFS)
    elseif iST == 3
        title('Example Stimuli and Spike Rasters (Background OFF)','FontSize',TitleFS)
    end
    
end


% Plot the PID through Time Results

% Set the offset for the different networks
OffSet = 10;
OffSet = [-OffSet,OffSet];

yMax = -inf;
yMin = inf;
for iST = 1:2
    for iPID = 1:4
        yMax = max([max(max(DP3.meanerror{iST}(:,:,iPID))),yMax]);
        yMin = min([min(min(DP3.meanerror{iST}(:,:,iPID))),yMin]);
    end
end

dy = yMax - yMin;
yMax = yMax + 0.1*dy;
yMin = yMin - 0.1*dy;

xMin = min([DP3.time{1}(1) - 0.5*(DP3.time{1}(2) - DP3.time{1}(1)),DP3.time{2}(1) - 0.5*(DP3.time{2}(2) - DP3.time{2}(1))]);
xMax = max([DP3.time{1}(end) - 0.5*(DP3.time{1}(2) - DP3.time{1}(1)),DP3.time{2}(end) - 0.5*(DP3.time{2}(2) - DP3.time{2}(1))]);

for iPID = 1:4

    f = subplot('Position',[LeftCoord(iPID),BottomCoord(3),width(iPID),height(3)]);
    hold on
    
    line([0,0],[yMin,yMax],'Color',PColor{3},'LineStyle','--','LineWidth',StimLineLW)
    text(-15,mean([yMin,yMax]),'Stim On','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',90,'FontSize',StimLabelFS,'Color',PColor{3})
    line([500,500],[yMin,yMax],'Color',PColor{2},'LineStyle','--','LineWidth',StimLineLW)
    text(500 + 15,mean([yMin,yMax]),'Stim Off','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-90,'FontSize',StimLabelFS,'Color',PColor{2})
    
    
    for iST = 1:2
        nBins = size(DP3.meanerror{iST},2);
        
        for iBin = 1:nBins
            line([DP3.time{iST}(iBin) + OffSet(iST),DP3.time{iST}(iBin) + OffSet(iST)],[DP3.meanerror{iST}(4,iBin,iPID),DP3.meanerror{iST}(5,iBin,iPID)],...
                'Color',PColor{iST + 3},'LineWidth',BPThinLW)
            scatter(DP3.time{iST}(iBin) + OffSet(iST),DP3.meanerror{iST}(1,iBin,iPID),BPDotSz,PColor{iST + 3},'filled')
            
            %         line([DP2.time{iST}(iBin),DP2.time{iST}(iBin)],[DP2.boxpoints{iST}(1,iBin),DP2.boxpoints{iST}(5,iBin)],'Color','k','LineWidth',BPThinLW)
            %         line([DP2.time{iST}(iBin),DP2.time{iST}(iBin)],[DP2.boxpoints{iST}(2,iBin),DP2.boxpoints{iST}(4,iBin)],'Color','k','LineWidth',BPThickLW)
            %         scatter(DP2.time{iST}(iBin),DP2.boxpoints{iST}(3,iBin),BPDotSz,PColor{2},'d','filled')
        end
    end
    
    ylim([yMin,yMax])
    xlim([xMin,xMax])
    
    
    set(gca,'FontSize',UnitFS)
    xlabel('Time (ms)','FontSize',AxLabelFS)
    ylabel([PIDNamesShort{iPID},' (bits)'],'FontSize',AxLabelFS)
    title(PIDNamesLong{iPID},'FontSize',TitleFS)
    
    if iPID == 1
        x1 = xMin + 0.4*(xMax - xMin);
        x2 = xMin + 0.43*(xMax - xMin);
        y1 = yMin + 0.9*(yMax - yMin);
        y2 = yMin + 0.8*(yMax - yMin);
        scatter(x1,y1,BPDotSz,PColor{4},'filled')
        text(x2,y1,'BG Ex ON','FontSize',LegFS)
        scatter(x1,y2,BPDotSz,PColor{5},'filled')
        text(x2,y2,'BG Ex OFF','FontSize',LegFS)
    end
    
end


% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo4Nor')
cd(SimDir)


%% Make the figure for the XOR gate

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.18;

% Set the top and bottom margins in inches
tmargin = 0.23;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace = [0.4,0.4,0.4];

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set Paper size
papersize = [7.4 3.5];

% Set the figure panel dimensions (number of horizontal figures by number
% of vertical figures)
figdim = [4,2];

% Set the height proportions for each row
heightprop = [1,1];

% Set the width proportions for each column
widthprop = [1,1,1,1];

% calculate the width and height based on dimensions above in inches
widthtot = papersize(1) - lmargin - rmargin - sum(hspace);
widthunit = widthtot/sum(widthprop);
width = widthunit*widthprop;
heighttot = papersize(2) - tmargin - bmargin - (figdim(2) - 1)*vspace;
heightunit = heighttot/sum(heightprop);
height = heightunit*heightprop;

LeftCoord = cumsum([0,width(1:(end - 1)) + hspace]) + lmargin;
BottomCoord = cumsum([0,fliplr(height + vspace)]) + bmargin;
BottomCoord(end) = [];
BottomCoord = fliplr(BottomCoord);

% Convert to fraction of page sizes
LeftCoord = LeftCoord/papersize(1);
BottomCoord = BottomCoord/papersize(2);
width = width/papersize(1);
height = height/papersize(2);
hspace = hspace/papersize(1);

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
StimLabelFS = 6;

% Set the linewidths
VoltageLW = 0.6;
BPThickLW = 3;
BPThinLW = 0.6;
StimLineLW = 1;
DelayLW = 2;
AxLW = 0.4;
SpikeLW = 0.4;

% Set dot sizes
BPDotSz = 8;

% Make PID names
PIDNamesLong = {'Redundancy','Unique A','Unique B','Synergy'};
PIDNamesShort = {'Red','Uniq A','Uniq B','Syn'};

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');

% Make the spike rasters

iST = 4;

f = subplot('Position',[LeftCoord(3),BottomCoord(1),width(3) + width(4) + hspace(3),height(1)]);
hold on

% Set the limits
xMin = 0;
xMax = max(DP1.time{iST});
yMax = max([DP1.current{iST,1},DP1.current{iST,2}]);
yMin = min([DP1.current{iST,1},DP1.current{iST,2}]);
dy = yMax - yMin;

% Plot the stimulus current
plot(DP1.time{iST},DP1.current{iST,1},'k','LineWidth',VoltageLW)
plot(DP1.time{iST},DP1.current{iST,2} + yMax + dy,'k','LineWidth',VoltageLW) % Assume inputs have same magnitude

text(xMin + 0.18*(xMax - xMin),0.5*dy,'Stimulus A','FontSize',StimLabelFS)
text(xMin + 0.18*(xMax - xMin),2.5*dy,'Stimulus B','FontSize',StimLabelFS)

line([xMin + 0.1*(xMax - xMin),xMin + 0.1*(xMax - xMin)],[0.5*dy,1.5*dy],'Color','k','LineWidth',VoltageLW)
text(xMin + 0.09*(xMax - xMin),dy,[num2str(dy),' pA'],'FontSize',LegFS,'HorizontalAlignment','right')

yMax = yMax + dy + yMax;

% Plot the spikes
for iSpike = 1:length(DP1.spkTimes{iST,1})
    line([DP1.spkTimes{iST,1}(iSpike),DP1.spkTimes{iST,1}(iSpike)],[yMax + dy,yMax + 2*dy],'Color',PColor{2},'LineWidth',SpikeLW)
end
for iSpike = 1:length(DP1.spkTimes{iST,2})
    line([DP1.spkTimes{iST,2}(iSpike),DP1.spkTimes{iST,2}(iSpike)],[yMax + 3*dy,yMax + 4*dy],'Color',PColor{3},'LineWidth',SpikeLW)
end
for iSpike = 1:length(DP1.spkTimes{iST,3})
    line([DP1.spkTimes{iST,3}(iSpike),DP1.spkTimes{iST,3}(iSpike)],[yMax + 5*dy,yMax + 6*dy],'Color',PColor{4},'LineWidth',SpikeLW)
end
for iSpike = 1:length(DP1.spkTimes{iST,4})
    line([DP1.spkTimes{iST,4}(iSpike),DP1.spkTimes{iST,4}(iSpike)],[yMax + 7*dy,yMax + 8*dy],'Color',PColor{1},'LineWidth',SpikeLW)
end

text(xMin + 0.18*(xMax - xMin),yMax + 2.5*dy,'Neuron E1 Spikes','FontSize',StimLabelFS,'Color',PColor{2})
text(xMin + 0.18*(xMax - xMin),yMax + 4.5*dy,'Neuron E2 Spikes','FontSize',StimLabelFS,'Color',PColor{3})
text(xMin + 0.18*(xMax - xMin),yMax + 6.5*dy,'Neuron I1 Spikes','FontSize',StimLabelFS,'Color',PColor{4})
text(xMin + 0.18*(xMax - xMin),yMax + 8.5*dy,'Neuron E3 Spikes','FontSize',StimLabelFS,'Color',PColor{1})

xlim([xMin,xMax])
ylim([yMin - dy,yMax + 9*dy])
set(gca,'FontSize',UnitFS)
xlabel('Time (ms)','FontSize',AxLabelFS)
title('Example Stimuli and Spike Rasters','FontSize',TitleFS)
set(gca,'ytick',[])


% Plot the PID through Time Results

yMax = -inf;
yMin = inf;
for iPID = 1:4
    yMax = max([max(max(DP4.meanerror(:,:,iPID))),yMax]);
    yMin = min([min(min(DP4.meanerror(:,:,iPID))),yMin]);
end

dy = yMax - yMin;
yMax = yMax + 0.1*dy;
yMin = yMin - 0.1*dy;

xMin = DP4.time(1) - 0.5*(DP4.time(2) - DP4.time(1));
xMax = DP4.time(end) + 0.5*(DP4.time(2) - DP4.time(1));

for iPID = 1:4

    f = subplot('Position',[LeftCoord(iPID),BottomCoord(2),width(iPID),height(2)]);
    hold on
    
    nBins = size(DP4.meanerror,2);
    
    line([0,0],[yMin,yMax],'Color',PColor{3},'LineStyle','--','LineWidth',StimLineLW)
    text(-15,mean([yMin,yMax]),'Stim On','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',90,'FontSize',StimLabelFS,'Color',PColor{3})
    line([500,500],[yMin,yMax],'Color',PColor{2},'LineStyle','--','LineWidth',StimLineLW)
    text(500 + 15,mean([yMin,yMax]),'Stim Off','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-90,'FontSize',StimLabelFS,'Color',PColor{2})
    
    for iBin = 1:nBins
        line([DP4.time(iBin),DP4.time(iBin)],[DP4.meanerror(4,iBin,iPID),DP4.meanerror(5,iBin,iPID)],'Color','k','LineWidth',BPThinLW)
        scatter(DP4.time(iBin),DP4.meanerror(1,iBin,iPID),BPDotSz,'k','filled')
    end
    
    ylim([yMin,yMax])
    xlim([xMin,xMax])
    
    
    set(gca,'FontSize',UnitFS)
    xlabel('Time (ms)','FontSize',AxLabelFS)
    ylabel([PIDNamesShort{iPID},' (bits)'],'FontSize',AxLabelFS)
    title(PIDNamesLong{iPID},'FontSize',TitleFS)
    
end



% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo4XOR')
cd(SimDir)


%% Make the figure for the OR gate

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.18;

% Set the top and bottom margins in inches
tmargin = 0.23;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace = [0.4,0.4,0.4];

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set Paper size
papersize = [7.4 4.75];

% Set the figure panel dimensions (number of horizontal figures by number
% of vertical figures)
figdim = [4,3];

% Set the height proportions for each row
heightprop = [1,1,1];

% Set the width proportions for each column
widthprop = [1,1,1,1];

% calculate the width and height based on dimensions above in inches
widthtot = papersize(1) - lmargin - rmargin - sum(hspace);
widthunit = widthtot/sum(widthprop);
width = widthunit*widthprop;
heighttot = papersize(2) - tmargin - bmargin - (figdim(2) - 1)*vspace;
heightunit = heighttot/sum(heightprop);
height = heightunit*heightprop;

LeftCoord = cumsum([0,width(1:(end - 1)) + hspace]) + lmargin;
BottomCoord = cumsum([0,fliplr(height + vspace)]) + bmargin;
BottomCoord(end) = [];
BottomCoord = fliplr(BottomCoord);

% Convert to fraction of page sizes
LeftCoord = LeftCoord/papersize(1);
BottomCoord = BottomCoord/papersize(2);
width = width/papersize(1);
height = height/papersize(2);
hspace = hspace/papersize(1);

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
StimLabelFS = 6;

% Set the linewidths
VoltageLW = 0.6;
BPThickLW = 3;
BPThinLW = 0.6;
StimLineLW = 1;
DelayLW = 2;
AxLW = 0.4;
SpikeLW = 0.4;

% Set dot sizes
BPDotSz = 8;

% Make PID names
PIDNamesLong = {'Redundancy','Unique A','Unique B','Synergy'};
PIDNamesShort = {'Red','Uniq A','Uniq B','Syn'};

% Make the correlation names
CorrNames = {'Anti-Correlated (a = -0.25)','Uncorrelated (a = 0)','Correlated (a = 0.25)'};

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');


% Make the spike rasters

iST = 5;

f = subplot('Position',[LeftCoord(3),BottomCoord(1),width(3) + width(4) + hspace(3),height(1)]);
hold on

% Set the limits
xMin = 0;
xMax = max(DP1.time{iST});
yMax = max([DP1.current{iST,1},DP1.current{iST,2}]);
yMin = min([DP1.current{iST,1},DP1.current{iST,2}]);
dy = yMax - yMin;

% Plot the stimulus current
plot(DP1.time{iST},DP1.current{iST,1},'k','LineWidth',VoltageLW)
plot(DP1.time{iST},DP1.current{iST,2} + yMax + dy,'k','LineWidth',VoltageLW) % Assume inputs have same magnitude

text(xMin + 0.18*(xMax - xMin),0.5*dy,'Stimulus A','FontSize',StimLabelFS)
text(xMin + 0.18*(xMax - xMin),2.5*dy,'Stimulus B','FontSize',StimLabelFS)

line([xMin + 0.1*(xMax - xMin),xMin + 0.1*(xMax - xMin)],[0.5*dy,1.5*dy],'Color','k','LineWidth',VoltageLW)
text(xMin + 0.09*(xMax - xMin),dy,[num2str(dy),' pA'],'FontSize',LegFS,'HorizontalAlignment','right')

yMax = yMax + dy + yMax;

% Plot the spikes
for iSpike = 1:length(DP1.spkTimes{iST,1})
    line([DP1.spkTimes{iST,1}(iSpike),DP1.spkTimes{iST,1}(iSpike)],[yMax + dy,yMax + 2*dy],'Color',PColor{2},'LineWidth',SpikeLW)
end
for iSpike = 1:length(DP1.spkTimes{iST,2})
    line([DP1.spkTimes{iST,2}(iSpike),DP1.spkTimes{iST,2}(iSpike)],[yMax + 3*dy,yMax + 4*dy],'Color',PColor{3},'LineWidth',SpikeLW)
end
for iSpike = 1:length(DP1.spkTimes{iST,4})
    line([DP1.spkTimes{iST,4}(iSpike),DP1.spkTimes{iST,4}(iSpike)],[yMax + 5*dy,yMax + 6*dy],'Color',PColor{1},'LineWidth',SpikeLW)
end

text(xMin + 0.18*(xMax - xMin),yMax + 2.5*dy,'Neuron E1 Spikes','FontSize',StimLabelFS,'Color',PColor{2})
text(xMin + 0.18*(xMax - xMin),yMax + 4.5*dy,'Neuron E2 Spikes','FontSize',StimLabelFS,'Color',PColor{3})
text(xMin + 0.18*(xMax - xMin),yMax + 6.5*dy,'Neuron E3 Spikes','FontSize',StimLabelFS,'Color',PColor{1})

xlim([xMin,xMax])
ylim([yMin - dy,yMax + 7*dy])
set(gca,'FontSize',UnitFS)
xlabel('Time (ms)','FontSize',AxLabelFS)
title('Example Stimuli and Spike Rasters','FontSize',TitleFS)
set(gca,'ytick',[])


% Plot the PID through Time Results

yMax = -inf;
yMin = inf;
for iPID = [1,4]
    for iST = 1:3
        yMax = max([max(max(DP5.meanerror{iST}(:,:,iPID))),yMax]);
        yMin = min([min(min(DP5.meanerror{iST}(:,:,iPID))),yMin]);
    end
end

dy = yMax - yMin;
yMax = yMax + 0.1*dy;
yMin = yMin - 0.1*dy;

xMin = min([DP5.time{1}(1) - 0.5*(DP5.time{1}(2) - DP5.time{1}(1)),...
    DP5.time{2}(1) - 0.5*(DP5.time{2}(2) - DP5.time{2}(1)),DP5.time{3}(1) - 0.5*(DP5.time{3}(2) - DP5.time{3}(1))]);
xMax = max([DP5.time{1}(end) + 0.5*(DP5.time{1}(2) - DP5.time{1}(1)),...
    DP5.time{2}(end) + 0.5*(DP5.time{2}(2) - DP5.time{2}(1)),DP5.time{3}(end) + 0.5*(DP5.time{3}(2) - DP5.time{3}(1))]);

for iPID = [1,4]
    for iST = 1:3
        
        if iPID == 1
            iRow = 2;
        elseif iPID == 4
            iRow = 3;
        end
        
        f = subplot('Position',[LeftCoord(iST),BottomCoord(iRow),width(iST),height(iRow)]);
        hold on
        
        nBins = size(DP5.meanerror{iST},2);
        
        line([0,0],[yMin,yMax],'Color',PColor{3},'LineStyle','--','LineWidth',StimLineLW)
        text(-15,mean([yMin,yMax]),'Stim On','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',90,'FontSize',StimLabelFS,'Color',PColor{3})
        line([500,500],[yMin,yMax],'Color',PColor{2},'LineStyle','--','LineWidth',StimLineLW)
        text(500 + 15,mean([yMin,yMax]),'Stim Off','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-90,'FontSize',StimLabelFS,'Color',PColor{2})
        
        for iBin = 1:nBins
            line([DP5.time{iST}(iBin),DP5.time{iST}(iBin)],[DP5.meanerror{iST}(4,iBin,iPID),DP5.meanerror{iST}(5,iBin,iPID)],'Color','k','LineWidth',BPThinLW)
            scatter(DP5.time{iST}(iBin),DP5.meanerror{iST}(1,iBin,iPID),BPDotSz,'k','filled')
        end
        
        ylim([yMin,yMax])
        xlim([xMin,xMax])
        
        
        set(gca,'FontSize',UnitFS)
        xlabel('Time (ms)','FontSize',AxLabelFS)
        ylabel([PIDNamesShort{iPID},' (bits)'],'FontSize',AxLabelFS)
        title([PIDNamesLong{iPID},char(10),CorrNames{iST}],'FontSize',TitleFS)
        
    end
end


% Plot the PID as a function of correlation results

yMax = -inf;
yMin = inf;
for iPID = [1,4]
    yMax = max([max(max(DP6.InhibVar(6:10,:,iPID))),yMax]);
    yMin = min([min(min(DP6.InhibVar(6:10,:,iPID))),yMin]);
end

dy = yMax - yMin;
yMax = yMax + 0.1*dy;
yMin = yMin - 0.1*dy;

xMin = DP6.CorParam(1);
xMax = DP6.CorParam(end);

for iPID = [1,4]
    
    if iPID == 1
        iRow = 2;
    elseif iPID == 4
        iRow = 3;
    end
    
    f = subplot('Position',[LeftCoord(4),BottomCoord(iRow),width(4),height(iRow)]);
    hold on
    
    nBins = length(DP6.CorParam);
    
    for iBin = 1:nBins
        line([DP6.CorParam(iBin),DP6.CorParam(iBin)],[DP6.InhibVar(9,iBin,iPID),DP6.InhibVar(10,iBin,iPID)],'Color','k','LineWidth',BPThinLW)
        scatter(DP6.CorParam(iBin),DP6.InhibVar(6,iBin,iPID),BPDotSz,'k','filled')
    end
    
    ylim([yMin,yMax])
    xlim([xMin,xMax])
    
    
    set(gca,'FontSize',UnitFS)
    xlabel('a','FontSize',AxLabelFS)
    ylabel([PIDNamesShort{iPID},' (bits)'],'FontSize',AxLabelFS)
    title([PIDNamesLong{iPID},' Input',char(10),'Correlation Dependence'],'FontSize',TitleFS)
    
end


% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo4OR')
cd(SimDir)