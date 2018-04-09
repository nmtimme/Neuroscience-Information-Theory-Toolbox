%% neurDemo2
%
%   neuroDemo2 reproduces the 2 or 3 neuron information transmition
%   figures. Note, the scratch directory in the simulations folder will be 
%   used to save data. These simulations will take up about 2 minutes to 
%   run and it will require about 12 GB of space.
%
% Other m-files required: modeltype2, instinfo, inverseFNoise,
% simpleModels_FSI, simpleModel_RS
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% October 2016; Last revision: 28-Oct-2016

%% Set Parameters

disp('Note that this simulation will take about 15 minutes to run and it requires about 12 GB of space.')

% Set the number of models to run for each subtype
nModelsPerSubtype = 20;

% Set the number of weights to scan for the changing encoding
nScan = 11;

% Note, most parameters for the simulations are set in modeltype2

% Set the weights for the model subtypes
weights = 200*ones([nScan + 1,3]);
weights(1:nScan,3) = linspace(0,-150,nScan)';
weights(nScan + 1,3) = 0;

% Set the noiseScales for the model subtypes
noiseScale = [0.2*ones([1,nScan]),0.1];

% Find the number of subtypes
nSTs = nScan + 1;

% Set the pulseSize
pulseSize = [500,0];

% Set the bin size in milliseconds for the different subtypes
BinSizeA = 50; % All but last subtype
BinSizeB = 25; % Last subtype

% Set the time window around the start of the stimuli to analyze in
% milliseconds
MaxLead = 200;
MaxLag = 800;


% Get the scratch directory to use to store data
if ispc
    ScratchDir = [pwd,'\Scratch\ModelType2\'];
elseif isunix
    ScratchDir = [pwd,'/Scratch/ModelType2/'];
elseif ismac
    ScratchDir = [pwd,'/Scratch/ModelType2/'];
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
        
        [mV,spkTimes,memNoise,current,stim,modelParams] = modeltype2(pulseSize, weights(iSubType,:), 'noiseScale', noiseScale(iSubType)*ones([1,3]));
        
        save([ScratchDir,'ModelType2ST',num2str(iSubType),'Run',num2str(iModel),'.mat'],'mV','spkTimes','memNoise','stim','modelParams')
        
        disp(['Finished Generating Data for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Convert the Raw Data to DataRaster Format

% Go to the directory with the analysis software
cd(AnaDir)

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType2ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])
        
        % Find the number of neurons
        nNeurons = size(spkTimes,1);
        
        % Figure out which bin size to use
        if iSubType ~= nSTs
            BinSize = BinSizeA;
        else
            BinSize = BinSizeB;
        end
        
        % Make the Data Raster
        DataRaster = cell([2,1]);
        
        % Format the Data
        [Temp,timeboundaries] = formattool(ones(size(spkTimes{1})),spkTimes{1},stim.times,BinSize,MaxLead,MaxLag,'count');
        DataRaster{1} = NaN([nNeurons,size(Temp,2),size(Temp,3)]);
        DataRaster{1}(1,:,:) = Temp;
        for i = 2:nNeurons
            DataRaster{1}(i,:,:) = formattool(ones(size(spkTimes{i})),spkTimes{i},stim.times,BinSize,MaxLead,MaxLag,'count');
        end
        DataRaster{2} = reshape(stim.types,[1,1,length(stim.types)]);
        
        % Record the binsize
        modelParams.binsize = BinSize;
        modelParams.maxlead = MaxLead;
        modelParams.maxlag = MaxLag;
        modelParams.Time = timeboundaries;
        
        % Save the Data
        save([ScratchDir,'ModelType2ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'],'DataRaster','modelParams')
        
        disp(['Finished Converting Data to DataRaster Format for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Convert the DataRaster Format Data to StatesRasters

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType2ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])
        
        % Make the stating method assignment array (making sure to count the
        % number of pulses that were used)
        nPulses = length(modelParams.pulseSize);
        MethodAssign = {2,1,'Nat',{};...
            1,1,'UniCB',{nPulses};...
            1,2,'UniCB',{nPulses};...
            1,3,'UniCB',{nPulses}};
        
        % Convert the data
        StatesRaster = data2states(DataRaster,MethodAssign);
        
        % Save the Data
        save([ScratchDir,'ModelType2ST',num2str(iSubType),'Run',num2str(iModel),'SR.mat'],'StatesRaster','modelParams')
        
        disp(['Finished Converting Data to StatesRaster Format for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Perform the Information Theory Analysis

InfoResults = cell([nSTs,nModelsPerSubtype]);
SubTypeParams = cell([nSTs,1]);
for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType2ST',num2str(iSubType),'Run',num2str(iModel),'SR.mat'])
        
        % Save the subtype parameters. We assume all runs of a subtype used the
        % same parameters.
        if iModel == 1
            SubTypeParams{iSubType} = modelParams;
        end
        
        % Figure out how many time steps we have in this data set
        nT = size(StatesRaster{1},2);
        
        % Perform the appropriate analyses based on the subtype
        if iSubType ~= nSTs
            
            % Preallocate space
            MIE2Vals = NaN([1,nT]);
            MIE1Vals = NaN([1,nT]);
            MII1E2Vals = NaN([1,nT]);
            
            % Stimulus encoding by E2
            Method = 'PairMI';
            for iT = 1:nT
                VariableIDs = {1,3,iT;2,1,1};
                MIE2Vals(iT) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            % Stimulus encoding by E1
            Method = 'PairMI';
            for iT = 1:nT
                VariableIDs = {1,1,iT;2,1,1};
                MIE1Vals(iT) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            % Stimulus encoding by {I1,E2}
            Method = 'JointMI2';
            for iT = 1:nT
                VariableIDs = {1,2,iT;1,3,iT;2,1,1};
                MII1E2Vals(iT) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            % Put the results in the correct variable
            InfoResults{iSubType,iModel} = [MIE2Vals;MIE1Vals;MII1E2Vals];
            
        else
            
            % Preallocate space
            MIE2Vals = NaN([1,nT]);
            TEVals = NaN([1,nT]);
            ITVals = NaN([1,nT]);
            
            % Stimulus encoding by E2
            Method = 'PairMI';
            for iT = 1:nT
                VariableIDs = {1,3,iT;2,1,1};
                MIE2Vals(iT) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            % TE from E1 to E2
            Method = 'TE';
            for iT = 2:nT
                VariableIDs = {1,3,iT;1,3,iT - 1;1,1,iT - 1};
                TEVals(iT) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            % IT from E1 to E2
            Method = 'InfoTrans';
            for iT = 2:nT
                VariableIDs = {1,3,iT;1,3,iT - 1;2,1,1;1,1,iT - 1};
                ITVals(iT) = instinfo(StatesRaster, Method, VariableIDs);
            end
            
            % Put the results in the correct variable
            InfoResults{iSubType,iModel} = [MIE2Vals;TEVals;ITVals];
            
        end
        
        disp(['Finished Analyzing for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

% Save the Results
save([ScratchDir,'ModelType2InfoResults.mat'],'InfoResults','SubTypeParams')

%% Gather Data for Plotting

% Get Example Voltage Traces

DP1 = struct;

DP1.spkTimes = cell([nSTs,3]);
DP1.time = cell([nSTs,1]);
DP1.mV = cell([nSTs,3]);
DP1.current = cell([nSTs,1]);
for iST = 1:nSTs
    
    % Load the first run
    load([ScratchDir,'ModelType2ST',num2str(iST),'Run1.mat'])
    
    % Get the total time vector
    time = (1:length(mV))*modelParams.tau;
    
    % Set the minimum and maximum times based on the subtype
    stimset = 10;
    
    ExampleTime = length(modelParams.pulseSize)*modelParams.intrastimT + (length(modelParams.pulseSize) + 1)*modelParams.interstimT;
    OneSetTime = length(modelParams.pulseSize)*modelParams.intrastimT + length(modelParams.pulseSize)*modelParams.interstimT;
    MinTime = (stimset - 1)*OneSetTime;
    MaxTime = MinTime + ExampleTime;
    
    % Get the time vector
    DP1.time{iST} = time((time > MinTime) & (time <= MaxTime));
    
    for i = 1:3
        
        % Get the spike times,
        DP1.spkTimes{iST,i} = spkTimes{i}(ismember(spkTimes{i},DP1.time{iST}));
        
        % Get the membrane voltage
        DP1.mV{iST,i} = mV(i,(time > MinTime) & (time <= MaxTime));
        
    end
    
    % Get the stimulus current
    DP1.current{iST} = stim.currentOriginal((time > MinTime) & (time <= MaxTime));
    
    % Correct offset in time for plotting
    DP1.time{iST} = DP1.time{iST} - MinTime;
    for i = 1:3
        DP1.spkTimes{iST,i} = DP1.spkTimes{iST,i} - MinTime;
    end
end

% Get the information results from the no inhibition model

DP2 = struct;

iSubType = nSTs;

modelParams = SubTypeParams{iSubType};

% Get the time values
DP2.time = mean(modelParams.Time,2)';
nBins = length(DP2.time);

% Get all the information values
AllInfo = zeros([nModelsPerSubtype,nBins,3]);
for iModel = 1:nModelsPerSubtype
    AllInfo(iModel,:,1) = InfoResults{iSubType,iModel}(1,:);
    AllInfo(iModel,:,2) = InfoResults{iSubType,iModel}(2,:);
    AllInfo(iModel,:,3) = InfoResults{iSubType,iModel}(3,:);
end

DP2.boxpoints = NaN([5,nBins,3]);
DP2.meanerror = NaN([5,nBins,3]);
for iBin = 1:nBins
    for iMeas = 1:3
        temp = sort(AllInfo(:,iBin,iMeas),'descend');
        DP2.boxpoints(1,iBin,iMeas) = temp(1);
        DP2.boxpoints(2,iBin,iMeas) = temp(floor(nModelsPerSubtype*0.75));
        DP2.boxpoints(3,iBin,iMeas) = median(temp);
        DP2.boxpoints(4,iBin,iMeas) = temp(ceil(nModelsPerSubtype*0.25));
        DP2.boxpoints(5,iBin,iMeas) = temp(end);
        DP2.meanerror(1,iBin,iMeas) = mean(temp);
        DP2.meanerror(2,iBin,iMeas) = mean(temp) + std(temp)/sqrt(length(temp));
        DP2.meanerror(3,iBin,iMeas) = mean(temp) - std(temp)/sqrt(length(temp));
        DP2.meanerror(4,iBin,iMeas) = mean(temp) + std(temp);
        DP2.meanerror(5,iBin,iMeas) = mean(temp) - std(temp);
    end
end

% Get the information results for the variable inhibition model

DP3 = struct;

DP3.time = cell([nSTs - 1,1]);
DP3.boxpoints = cell([nSTs - 1,1]);
DP3.nsig = cell([nSTs - 1,1]);
DP3.meanerror = cell([nSTs - 1,1]);
DP3.InhibVar = NaN([10,nSTs - 1,3]);
DP3.InhibVarI = NaN([1,nSTs - 1]);
for iSubType = 1:(nSTs - 1)
    
    modelParams = SubTypeParams{iSubType};
    
    
    % Get the time values
    DP3.time{iSubType} = mean(modelParams.Time,2)';
    nBins = length(DP3.time{iSubType});
    
    % Get all the information values and p-values
    AllInfo = zeros([nModelsPerSubtype,nBins,3]);
    AllpValues = NaN([nModelsPerSubtype,nBins,3]);
    for iModel = 1:nModelsPerSubtype
        AllInfo(iModel,:,1) = InfoResults{iSubType,iModel}(1,:);
        AllInfo(iModel,:,2) = InfoResults{iSubType,iModel}(2,:);
        AllInfo(iModel,:,3) = InfoResults{iSubType,iModel}(3,:);
    end
    
    % Get the data for information as a function of time
    DP3.boxpoints{iSubType} = NaN([5,nBins,3]);
    DP3.meanerror{iSubType} = NaN([5,nBins,3]);
    for iBin = 1:nBins
        for iMeas = 1:3
            temp = sort(AllInfo(:,iBin,iMeas),'descend');
            DP3.boxpoints{iSubType}(1,iBin,iMeas) = temp(1);
            DP3.boxpoints{iSubType}(2,iBin,iMeas) = temp(floor(nModelsPerSubtype*0.75));
            DP3.boxpoints{iSubType}(3,iBin,iMeas) = median(temp);
            DP3.boxpoints{iSubType}(4,iBin,iMeas) = temp(ceil(nModelsPerSubtype*0.25));
            DP3.boxpoints{iSubType}(5,iBin,iMeas) = temp(end);
            DP3.meanerror{iSubType}(1,iBin,iMeas) = mean(temp);
            DP3.meanerror{iSubType}(2,iBin,iMeas) = mean(temp) + std(temp)/sqrt(length(temp));
            DP3.meanerror{iSubType}(3,iBin,iMeas) = mean(temp) - std(temp)/sqrt(length(temp));
            DP3.meanerror{iSubType}(4,iBin,iMeas) = mean(temp) + std(temp);
            DP3.meanerror{iSubType}(5,iBin,iMeas) = mean(temp) - std(temp);
        end
    end
    
    
    % Get average information during stimulus as a function of inhibition
    for iMeas = 1:3
        temp = AllInfo(:,(DP3.time{iSubType} > 0) & (DP3.time{iSubType} < modelParams.intrastimT),iMeas);
        temp = sort(temp(:),'descend');
        DP3.InhibVar(1,iSubType,iMeas) = temp(1);
        DP3.InhibVar(2,iSubType,iMeas) = temp(floor(nModelsPerSubtype*0.75));
        DP3.InhibVar(3,iSubType,iMeas) = median(temp);
        DP3.InhibVar(4,iSubType,iMeas) = temp(ceil(nModelsPerSubtype*0.25));
        DP3.InhibVar(5,iSubType,iMeas) = temp(end);
        DP3.InhibVar(6,iSubType,iMeas) = mean(temp);
        DP3.InhibVar(7,iSubType,iMeas) = mean(temp) + std(temp)/sqrt(length(temp));
        DP3.InhibVar(8,iSubType,iMeas) = mean(temp) - std(temp)/sqrt(length(temp));
        DP3.InhibVar(9,iSubType,iMeas) = mean(temp) + std(temp);
        DP3.InhibVar(10,iSubType,iMeas) = mean(temp) - std(temp);
    end
    
    % Record the inhibitory strength
    DP3.InhibVarI(iSubType) = modelParams.weights(3);
    
    
    
end

% Save the plotting results
save([ScratchDir,'ModelType2PlotData'],'DP1','DP2','DP3')

disp('Finished Organizing the Data for Plotting')

%% Make the two neuron figure

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.15;

% Set the top and bottom margins in inches
tmargin = 0.23;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace = [0.6,0.6];

% Set the vertical distance between plots in inches
vspace = 0.5;

% Set Paper size
papersize = [7.4 4];

% Set the figure panel dimensions (number of horizontal figures by number
% of vertical figures)
figdim = [3,2];

% Set the height proportions for each row
heightprop = [1,1];

% Set the width proportions for each column
widthprop = [1,1,1];

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

% Set Some Colors
PColor = {[0,0,1],[1,0,0],[0,210/255,50/255],[200/255,0,255/255],[238/255,130/255,238/255],[200/255,0,255/255],[0/255,255/255,185/255],[238/255,130/255,238/255]};
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
TinyLabelFS = 6;
LabelFS = 16;
InsetFS = 4;
ExpTextFS = 10;
DelayFS = 6;
StimLabelFS = 6;
DiaFS = 4;

% Set the linewidths
VoltageLW = 0.6;
BPThickLW = 3;
BPThinLW = 1;
StimLineLW = 1;
DelayLW = 2;
AxLW = 0.4;

% Set dot sizes
BPDotSz = 10;

% Set several quantities related to the second y-axis ticks
cTickL = 0.03;
cTickspacing = 0.03;
cTickLabelspacing = [0.22,0.25]; 

% Set parameters for the diagrams
NeuriteLW = 0.7;
NeurDotSz = 100;
NeurXLoc = 0.6;
yPC1 = 0.03;
yPC2 = 1;
PCDeg = 30;

% Set the model subtype
iST = length(DP1.current);


% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');


% Make the spike rasters

f = subplot('Position',[LeftCoord(2),BottomCoord(1),width(2) + width(3) + hspace(2),height(1)]);
hold on

% Set the limits
xMin = 0;
xMax = max(DP1.time{iST});
yMax = max(DP1.current{iST});
yMin = min(DP1.current{iST});
dy = yMax - yMin;

% Plot the stimulus current
plot(DP1.time{iST},DP1.current{iST},'k','LineWidth',VoltageLW)

% Plot the spikes
for iSpike = 1:length(DP1.spkTimes{iST,1})
    line([DP1.spkTimes{iST,1}(iSpike),DP1.spkTimes{iST,1}(iSpike)],[yMax + 0.1*dy,yMax + 0.2*dy],'Color',PColor{2},'LineWidth',VoltageLW)
end
for iSpike = 1:length(DP1.spkTimes{iST,3})
    line([DP1.spkTimes{iST,3}(iSpike),DP1.spkTimes{iST,3}(iSpike)],[yMax + 0.35*dy,yMax + 0.45*dy],'Color',PColor{1},'LineWidth',VoltageLW)
end

text(xMin + 0.6*(xMax - xMin),yMax + 0.05*dy,'Neuron E1 Spikes','FontSize',StimLabelFS,'Color',PColor{2})
text(xMin + 0.6*(xMax - xMin),yMax + 0.275*dy,'Neuron E2 Spikes','FontSize',StimLabelFS,'Color',PColor{1})

xlim([xMin,xMax])
ylim([yMin - 0.1*dy,yMax + 0.5*dy])
set(gca,'FontSize',UnitFS)
xlabel('Time (ms)','FontSize',AxLabelFS)
ylabel('Stimulus Current (pA)','FontSize',AxLabelFS)
title('Example Stimulus Trace and Spike Raster','FontSize',TitleFS)






% Plot the Information

% Make uniform yAxis
yMax = max(max(DP2.meanerror(:)));
yMin = min(min(DP2.meanerror(:)));
dy = yMax - yMin;
yMax = yMax + 0.1*dy;
yMin = yMin - 0.1*dy;

for iMeas = 1:3
    
    f = subplot('Position',[LeftCoord(iMeas),BottomCoord(2),width(iMeas),height(2)]);
    hold on
    
    nBins = size(DP2.meanerror,2);
    
    line([0,0],[yMin,yMax],'Color',PColor{3},'LineStyle','--','LineWidth',StimLineLW)
    text(-15,mean([yMin,yMax]),'Stim On','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',90,'FontSize',StimLabelFS,'Color',PColor{3})
    line([500,500],[yMin,yMax],'Color',PColor{2},'LineStyle','--','LineWidth',StimLineLW)
    text(500 + 15,mean([yMin,yMax]),'Stim Off','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-90,'FontSize',StimLabelFS,'Color',PColor{2})
    
    for iBin = 1:nBins
        line([DP2.time(iBin),DP2.time(iBin)],[DP2.meanerror(4,iBin,iMeas),DP2.meanerror(5,iBin,iMeas)],'Color','k','LineWidth',BPThinLW)
        scatter(DP2.time(iBin),DP2.meanerror(1,iBin,iMeas),BPDotSz,'k','filled')
    end
    
    ylim([yMin,yMax])
    xlim([-200,800])
    
    
    set(gca,'FontSize',UnitFS)
    xlabel('Time (ms)','FontSize',AxLabelFS)
    if iMeas == 1
        ylabel('MI (bits)','FontSize',AxLabelFS)
        title('E2 Stimulus Encoding','FontSize',TitleFS)
    elseif iMeas == 2
        ylabel('TE (bits)','FontSize',AxLabelFS)
        title('E1 to E2 Transfer Entropy','FontSize',TitleFS)
    elseif iMeas == 3
        ylabel('IT (bits)','FontSize',AxLabelFS)
        title('E1 to E2 Stim Info Transmission','FontSize',TitleFS)
    end
    
end


% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo2TwoNeuronPanel')
cd(SimDir)


%% Make the Three figure

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.4;
rmargin = 0.18;

% Set the top and bottom margins in inches
tmargin = 0.23;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace = [0.6,0.6,0.6];

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set Paper size
papersize = [7.4 7];

% Set the figure panel dimensions (number of horizontal figures by number
% of vertical figures)
figdim = [4,4];

% Set the height proportions for each row
heightprop = [1,1,1,1];

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

% Set Some Colors
PColor = {[0,0,1],[1,0,0],[0,210/255,50/255],[150/255,0,200/255],[238/255,130/255,238/255],[200/255,0,255/255],[0/255,255/255,185/255],[238/255,130/255,238/255]};
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
TinyLabelFS = 6;
LabelFS = 16;
InsetFS = 4;
ExpTextFS = 10;
DelayFS = 6;
StimLabelFS = 6;
DiaFS = 4;

% Set the linewidths
VoltageLW = 0.6;
BPThickLW = 3;
BPThinLW = 1;
StimLineLW = 1;
DelayLW = 2;
AxLW = 0.4;

% Set dot sizes
BPDotSz = 10;

% Set several quantities related to the second y-axis ticks
cTickL = 0.03;
cTickspacing = 0.03;
cTickLabelspacing = [0.22,0.25]; 

% Set parameters for the diagrams
NeuriteLW = 0.7;
NeurDotSz = 100;
NeurXLoc = 0.6;
yPC1 = 0.03;
yPC2 = 1;
PCDeg = 30;

% Set the model subtypes to use
nSTs = length(DP1.current);
STFocus = [1,6,nSTs - 1];


% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');

% Make the E2 mutual information vs. inhibition plot

iMeas = 1;

f = subplot('Position',[LeftCoord(3),BottomCoord(1),width(3),height(1)]);
hold on

x = abs(DP3.InhibVarI);

xMin = min(x);
xMax = max(x);
dx = xMax - xMin;
yMin = min(min(DP3.InhibVar(6:10,:,iMeas)));
yMax = max(max(DP3.InhibVar(6:10,:,iMeas)));
dy = yMax - yMin;

% Plot the data points
for iST = 1:length(x)
    line([x(iST),x(iST)],[DP3.InhibVar(9,iST,iMeas),DP3.InhibVar(10,iST,iMeas)],'Color',PColor{1},'LineWidth',BPThinLW)
    scatter(x(iST),DP3.InhibVar(6,iST,iMeas),BPDotSz,PColor{1},'filled')
end

xlim([xMin - 0.05*dx,xMax + 0.05*dx])
ylim([yMin - 0.1*dy,yMax + 0.1*dy])
set(gca,'FontSize',UnitFS)
xlabel('Maximum Inhibition Current (pA)','FontSize',AxLabelFS)
ylabel('Mutual Information (bits)','FontSize',AxLabelFS)
title('E2 Stimulus Encoding','FontSize',TitleFS)




% Make the E2 and I1 mutual information vs. inhibition plot

OffSet = 0.015;

f = subplot('Position',[LeftCoord(4),BottomCoord(1),width(4),height(1)]);
hold on

x = abs(DP3.InhibVarI);

xMin = min(x);
xMax = max(x);
dx = xMax - xMin;
yMin = min(min(min(DP3.InhibVar(6:10,:,[2,3]))));
yMax = max(max(max(DP3.InhibVar(6:10,:,[2,3]))));
dy = yMax - yMin;
yMin = yMin - 0.3*dy;
yMax = yMax + 0.1*dy;
dy = yMax - yMin;

% Plot the data points for the E1 encoding
iMeas = 2;
for iST = 1:length(x)
    line([x(iST) - OffSet*dx,x(iST) - OffSet*dx],[DP3.InhibVar(9,iST,iMeas),DP3.InhibVar(10,iST,iMeas)],'Color',PColor{2},'LineWidth',BPThinLW)
    scatter(x(iST) - OffSet*dx,DP3.InhibVar(6,iST,iMeas),BPDotSz,PColor{2},'filled')
end

% Plot the data points for the I1 and E2 encoding
iMeas = 3;
for iST = 1:length(x)
    line([x(iST) + OffSet*dx,x(iST) + OffSet*dx],[DP3.InhibVar(9,iST,iMeas),DP3.InhibVar(10,iST,iMeas)],'Color',PColor{4},'LineWidth',BPThinLW)
    scatter(x(iST) + OffSet*dx,DP3.InhibVar(6,iST,iMeas),BPDotSz,PColor{4},'filled')
end

x1 = xMin + 0.6*dx;
x2 = xMin + 0.65*dx;
y1 = yMin + 0.15*dy;
y2 = yMin + 0.07*dy;

scatter(x1,y1,BPDotSz,PColor{2},'filled')
text(x2,y1,'E1','FontSize',LegFS)

scatter(x1,y2,BPDotSz,PColor{4},'filled')
text(x2,y2,'\{I1,E2\}','FontSize',LegFS)


xlim([xMin - 0.05*dx,xMax + 0.05*dx])
ylim([yMin,yMax])
set(gca,'FontSize',UnitFS)
xlabel('Maximum Inhibition Current (pA)','FontSize',AxLabelFS)
ylabel('Mutual Information (bits)','FontSize',AxLabelFS)
title('Stimulus Encoding Comparison','FontSize',TitleFS)


% Make the spike rasters

for iST = 1:3
    
    f = subplot('Position',[LeftCoord(1),BottomCoord(iST + 1),width(1) + width(2) + hspace(1),height(iST + 1)]);
    hold on
    
    % Set the limits
    xMin = 0;
    xMax = max(DP1.time{STFocus(iST)});
    yMax = max(DP1.current{STFocus(iST)});
    yMin = min(DP1.current{STFocus(iST)});
    dy = yMax - yMin;
    
    % Plot the stimulus current
    plot(DP1.time{STFocus(iST)},DP1.current{STFocus(iST)},'k','LineWidth',VoltageLW)
    
    % Plot the spikes
    for iSpike = 1:length(DP1.spkTimes{STFocus(iST),1})
        line([DP1.spkTimes{STFocus(iST),1}(iSpike),DP1.spkTimes{STFocus(iST),1}(iSpike)],[yMax + 0.1*dy,yMax + 0.2*dy],'Color',PColor{2},'LineWidth',VoltageLW)
    end
    for iSpike = 1:length(DP1.spkTimes{STFocus(iST),2})
        line([DP1.spkTimes{STFocus(iST),2}(iSpike),DP1.spkTimes{STFocus(iST),2}(iSpike)],[yMax + 0.35*dy,yMax + 0.45*dy],'Color',PColor{3},'LineWidth',VoltageLW)
    end
    for iSpike = 1:length(DP1.spkTimes{STFocus(iST),3})
        line([DP1.spkTimes{STFocus(iST),3}(iSpike),DP1.spkTimes{STFocus(iST),3}(iSpike)],[yMax + 0.6*dy,yMax + 0.7*dy],'Color',PColor{1},'LineWidth',VoltageLW)
    end
    
    text(xMin + 0.6*(xMax - xMin),yMax + 0.05*dy,'Neuron E1 Spikes','FontSize',StimLabelFS,'Color',PColor{2})
    text(xMin + 0.6*(xMax - xMin),yMax + 0.275*dy,'Neuron I1 Spikes','FontSize',StimLabelFS,'Color',PColor{3})
    text(xMin + 0.6*(xMax - xMin),yMax + 0.525*dy,'Neuron E2 Spikes','FontSize',StimLabelFS,'Color',PColor{1})
    
    xlim([xMin,xMax])
    ylim([yMin - 0.1*dy,yMax + 0.75*dy])
    set(gca,'FontSize',UnitFS)
    xlabel('Time (ms)','FontSize',AxLabelFS)
    ylabel('Stimulus Current (pA)','FontSize',AxLabelFS)
    if iST == 1
        title(['Low Inhibition(I_{max} = ',num2str(abs(DP3.InhibVarI(STFocus(iST)))),' pA)'],'FontSize',TitleFS)
    elseif iST == 2
        title(['Medium Inhibition (I_{max} = ',num2str(abs(DP3.InhibVarI(STFocus(iST)))),' pA)'],'FontSize',TitleFS)
    elseif iST == 3
        title(['Strong Inhibition (I_{max} = ',num2str(abs(DP3.InhibVarI(STFocus(iST)))),' pA)'],'FontSize',TitleFS)
    end
    
end

% Plot the E2 encoding plots

iMeas = 1;

yMax = -inf;
yMin = inf;
for iST = 1:3
    yMax = max([max(max(DP3.meanerror{STFocus(iST)}(:,:,iMeas))),yMax]);
    yMin = min([min(min(DP3.meanerror{STFocus(iST)}(:,:,iMeas))),yMin]);
end

dy = yMax - yMin;
yMax = yMax + 0.1*dy;
yMin = yMin - 0.1*dy;

for iST = 1:3

    f = subplot('Position',[LeftCoord(3),BottomCoord(1 + iST),width(3),height(1 + iST)]);
    hold on
    
    nBins = size(DP3.meanerror{STFocus(iST)},2);
    
    line([0,0],[yMin,yMax],'Color',PColor{3},'LineStyle','--','LineWidth',StimLineLW)
    text(-15,mean([yMin,yMax]),'Stim On','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',90,'FontSize',StimLabelFS,'Color',PColor{3})
    line([500,500],[yMin,yMax],'Color',PColor{2},'LineStyle','--','LineWidth',StimLineLW)
    text(500 + 15,mean([yMin,yMax]),'Stim Off','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-90,'FontSize',StimLabelFS,'Color',PColor{2})
    
    for iBin = 1:nBins
        line([DP3.time{STFocus(iST)}(iBin),DP3.time{STFocus(iST)}(iBin)],[DP3.meanerror{STFocus(iST)}(4,iBin,iMeas),DP3.meanerror{STFocus(iST)}(5,iBin,iMeas)],'Color','k','LineWidth',BPThinLW)
        scatter(DP3.time{STFocus(iST)}(iBin),DP3.meanerror{STFocus(iST)}(1,iBin,iMeas),BPDotSz,'k','filled')
    end
    
    ylim([yMin,yMax])
    xlim([-200,800])
    
    
    set(gca,'FontSize',UnitFS)
    xlabel('Time (ms)','FontSize',AxLabelFS)
    ylabel('Mutual Information (bits)','FontSize',AxLabelFS)
    title('E2 Stimulus Encoding','FontSize',TitleFS)
    
end
    
    
    

% Plot the E1 and {I1,E2} encoding plots

OffSet = 0.01;

yMax = -inf;
yMin = inf;
for iST = 1:3
    yMax = max([max(max(max(DP3.meanerror{STFocus(iST)}(:,:,[2,3])))),yMax]);
    yMin = min([max(min(min(DP3.meanerror{STFocus(iST)}(:,:,[2,3])))),yMin]);
end

dy = yMax - yMin;
yMax = yMax + 0.1*dy;
yMin = yMin - 0.1*dy;

for iST = 1:3

    f = subplot('Position',[LeftCoord(4),BottomCoord(1 + iST),width(4),height(1 + iST)]);
    hold on
    
    nBins = size(DP3.meanerror{STFocus(iST)},2);
    
    xMin = min(DP3.time{STFocus(iST)});
    xMax = max(DP3.time{STFocus(iST)});
    dx = xMax - xMin;
    
    line([0,0],[yMin,yMax],'Color',PColor{3},'LineStyle','--','LineWidth',StimLineLW)
    text(-15,mean([yMin,yMax]),'Stim On','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',90,'FontSize',StimLabelFS,'Color',PColor{3})
    line([500,500],[yMin,yMax],'Color',PColor{2},'LineStyle','--','LineWidth',StimLineLW)
    text(500 + 15,mean([yMin,yMax]),'Stim Off','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-90,'FontSize',StimLabelFS,'Color',PColor{2})
    
    for iBin = 1:nBins
        line([DP3.time{STFocus(iST)}(iBin) - OffSet*dx,DP3.time{STFocus(iST)}(iBin) - OffSet*dx],[DP3.meanerror{STFocus(iST)}(4,iBin,2),DP3.meanerror{STFocus(iST)}(5,iBin,2)],'Color',PColor{2},'LineWidth',BPThinLW)
        scatter(DP3.time{STFocus(iST)}(iBin) - OffSet*dx,DP3.meanerror{STFocus(iST)}(1,iBin,2),BPDotSz,PColor{2},'filled')
        line([DP3.time{STFocus(iST)}(iBin) + OffSet*dx,DP3.time{STFocus(iST)}(iBin) + OffSet*dx],[DP3.meanerror{STFocus(iST)}(4,iBin,3),DP3.meanerror{STFocus(iST)}(5,iBin,3)],'Color',PColor{4},'LineWidth',BPThinLW)
        scatter(DP3.time{STFocus(iST)}(iBin) + OffSet*dx,DP3.meanerror{STFocus(iST)}(1,iBin,3),BPDotSz,PColor{4},'filled')
    end
    
    x1 = xMin + 0.8*dx;
    x2 = xMin + 0.84*dx;
    y1 = yMin + 1*dy;
    y2 = yMin + 0.91*dy;
    
    scatter(x1,y1,BPDotSz,PColor{2},'filled')
    text(x2,y1,'E1','FontSize',LegFS)
    
    scatter(x1,y2,BPDotSz,PColor{4},'filled')
    text(x2,y2,'\{I1,E2\}','FontSize',LegFS)
    
    ylim([yMin,yMax])
    xlim([xMin - 0.05*dx,xMax + 0.05*dx])
    
    
    set(gca,'FontSize',UnitFS)
    xlabel('Time (ms)','FontSize',AxLabelFS)
    ylabel('Mutual Information (bits)','FontSize',AxLabelFS)
    title('Stimulus Encoding Comparison','FontSize',TitleFS)
    
end



% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo2ThreeNeuronPanel')
cd(SimDir)



