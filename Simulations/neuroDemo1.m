%% neurDemo1
%
%   neuroDemo1 reproduces the single neuron stimulation encoding figure. 
%   Note, the scratch directory in the simulations folder will be used to 
%   save data. These simulations will take up about 2 minutes to run and it
%   will require about 2 GB of space.
%
% Other m-files required: modeltype1, instinfo, inverseFNoise,
% simpleModels_FSI, simpleModel_RS
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% October 2016; Last revision: 28-Oct-2016

%% Set Parameters

disp('Note that this simulation will take about 2 minutes to run and it requires about 2 GB of space.')

% Set the number of models to run for each subtype
nModelsPerSubtype = 20;

% Set the number of subtypes for this analysis (must be 4)
nSTs = 4;

% Note, most parameters for the simulations are set in modeltype1

% Set the pulseSizes for the four model subtypes
pulseSize = cell([1,nSTs]);
pulseSize{1} = 500;
pulseSize{2} = [500,250];
pulseSize{3} = 500;
pulseSize{4} = [100,500,800];

% Set the stimulation format labels
stimFormat = cell([1,nSTs]);
stimFormat{1} = 'ONvsOFF';
stimFormat{2} = 'StimType';
stimFormat{3} = 'Delay';
stimFormat{4} = 'NonLinear';

% Set the delay in milliseconds for the delay model subtype
delay = 150;

% Set the bin size in milliseconds
BinSize = 50;

% Set the time window around the start of the stimuli to analyze in
% milliseconds
MaxLead = 200;
MaxLag = 800;

% Set the p-value threshold
AnalysisParams = struct;
AnalysisParams.MCpThresh = 0.01;

% Set the maximum number of Monte Carlo trials to perform
AnalysisParams.MCnSamples = 1000;

% Get the scratch directory to use to store data
if ispc
    ScratchDir = [pwd,'\Scratch\ModelType1\'];
elseif isunix
    ScratchDir = [pwd,'/Scratch/ModelType1/'];
elseif ismac
    ScratchDir = [pwd,'/Scratch/ModelType1/'];
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
        
        if ismember(iSubType,[1,2,4])
            [mV,spkTimes,memNoise,stim,modelParams] = modeltype1(stimFormat{iSubType},pulseSize{iSubType});
        else
            [mV,spkTimes,memNoise,stim,modelParams] = modeltype1(stimFormat{iSubType},pulseSize{iSubType},delay);
        end
        
        save([ScratchDir,'ModelType1ST',num2str(iSubType),'Run',num2str(iModel),'.mat'],'mV','spkTimes','memNoise','stim','modelParams')
        
        disp(['Finished Generating Data for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Convert the Raw Data to DataRaster Format

% Go to the directory with the analysis software
cd(AnaDir)

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType1ST',num2str(iSubType),'Run',num2str(iModel),'.mat'])
        
        % Make the DataRaster
        DataRaster = cell([2,1]);
        
        % Format the Data
        [DataRaster{1},timeboundaries] = formattool(ones(size(spkTimes)),spkTimes,stim.times,BinSize,MaxLead,MaxLag,'count');
        DataRaster{2} = reshape(stim.types,[1,1,length(stim.types)]);
        
        % Record the binsize
        modelParams.binsize = BinSize;
        modelParams.maxlead = MaxLead;
        modelParams.maxlag = MaxLag;
        modelParams.Time = timeboundaries;
        
        % Save the Data
        save([ScratchDir,'ModelType1ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'],'DataRaster','modelParams')
        
        disp(['Finished Converting Data to DataRaster Format for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Convert the DataRaster Format Data to StatesRasters

for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType1ST',num2str(iSubType),'Run',num2str(iModel),'DR.mat'])
        
        % Make the stating method assignment array (making sure to count the
        % number of pulses that were used)
        nPulses = length(modelParams.pulseSize);
        MethodAssign = {2,1,'Nat',{};...
            1,1,'UniCB',{nPulses}};
        
        % Convert the data
        StatesRaster = data2states(DataRaster,MethodAssign);
        
        % Save the Data
        save([ScratchDir,'ModelType1ST',num2str(iSubType),'Run',num2str(iModel),'SR.mat'],'StatesRaster','modelParams')
        
        disp(['Finished Converting Data to StatesRaster Format for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

%% Perform the Information Theory Analysis

InfoResults = cell([nSTs,nModelsPerSubtype]);
SubTypeParams = cell([nSTs,1]);
for iSubType = 1:nSTs
    for iModel = 1:nModelsPerSubtype
        
        % Load the data
        load([ScratchDir,'ModelType1ST',num2str(iSubType),'Run',num2str(iModel),'SR.mat'])
        
        % Save the subtype parameters. We assume all runs of a subtype used the
        % same parameters.
        if iModel == 1
            SubTypeParams{iSubType} = modelParams;
        end
        
        % Figure out how many time steps we have in this data set and allocate
        % space
        nT = size(StatesRaster{1},2);
        MIVals = NaN([1,nT]);
        
        % Set Parameters
        Method = 'PairMI';
        
        % Loop through the time bins
        for iT = 1:nT
            VariableIDs = {1,1,iT;2,1,1};
            MIVals(iT) = instinfo(StatesRaster, Method, VariableIDs);
        end
        
        InfoResults{iSubType,iModel} = MIVals;
        
        disp(['Finished Analyzing for Subtype ',num2str(iSubType),' Run ', num2str(iModel)])
        
    end
end

% Save the Results
save([ScratchDir,'ModelType1InfoResults.mat'],'InfoResults','SubTypeParams')

%% Gather Data for Plotting

% Get Example Voltage Traces

DP1 = struct;

DP1.spkTimes = cell([nSTs,1]);
DP1.time = cell([nSTs,1]);
DP1.mV = cell([nSTs,1]);
DP1.current = cell([nSTs,1]);
DP1.delay = NaN([nSTs,1]);
for iST = 1:nSTs
    
    % Load the first run
    load([ScratchDir,'ModelType1ST',num2str(iST),'Run1.mat'])
    
    % Get the total time vector
    time = (1:length(mV))*modelParams.tau;
    
    % Set the minimum and maximum times based on the subtype
    if iST == 4
        stimset = 9;
    else
        stimset = 10;
    end
    
    ExampleTime = length(modelParams.pulseSize)*modelParams.intrastimT + (length(modelParams.pulseSize) + 1)*modelParams.interstimT;
    OneSetTime = length(modelParams.pulseSize)*modelParams.intrastimT + length(modelParams.pulseSize)*modelParams.interstimT;
    MinTime = (stimset - 1)*OneSetTime;
    MaxTime = MinTime + ExampleTime;
    
    % Get the time vector
    DP1.time{iST} = time((time > MinTime) & (time <= MaxTime));
    
    % Get the spike times
    DP1.spkTimes{iST} = spkTimes(ismember(spkTimes,DP1.time{iST}));
    
    % Get the membrane voltage
    DP1.mV{iST} = mV((time > MinTime) & (time <= MaxTime));
    
    % Get the stimulus current
    DP1.current{iST} = stim.currentOriginal((time > MinTime) & (time <= MaxTime));
    
    % Get the delay
    DP1.delay(iST) = modelParams.delay;
    
    % Correct offset in time for plotting
    DP1.time{iST} = DP1.time{iST} - MinTime;
    DP1.spkTimes{iST} = DP1.spkTimes{iST} - MinTime;
end



%% Get the mutual information results

DP2 = struct;

DP2.delay = NaN([nSTs,1]);
DP2.time = cell([nSTs,1]);
DP2.boxpoints = cell([nSTs,1]);
DP2.meanerror = cell([nSTs,1]);
for iSubType = 1:nSTs
    
    modelParams = SubTypeParams{iSubType};
    
    % Get the delay
    DP2.delay(iSubType) = modelParams.delay;
    
    % Get the time values
    DP2.time{iSubType} = mean(modelParams.Time,2)';
    nBins = length(DP2.time{iSubType});
    
    % Get all the information values
    AllMI = zeros([nModelsPerSubtype,nBins]);
    for iModel = 1:nModelsPerSubtype
        AllMI(iModel,:) = InfoResults{iSubType,iModel};
    end
    
    DP2.boxpoints{iSubType} = NaN([5,nBins]);
    DP2.meanerror{iSubType} = NaN([5,nBins]);
    for iBin = 1:nBins
        temp = sort(AllMI(:,iBin),'descend');
        DP2.boxpoints{iSubType}(1,iBin) = temp(1);
        DP2.boxpoints{iSubType}(2,iBin) = temp(floor(nModelsPerSubtype*0.75));
        DP2.boxpoints{iSubType}(3,iBin) = median(temp);
        DP2.boxpoints{iSubType}(4,iBin) = temp(ceil(nModelsPerSubtype*0.25));
        DP2.boxpoints{iSubType}(5,iBin) = temp(end);
        DP2.meanerror{iSubType}(1,iBin) = mean(temp);
        DP2.meanerror{iSubType}(2,iBin) = mean(temp) + std(temp)/sqrt(length(temp));
        DP2.meanerror{iSubType}(3,iBin) = mean(temp) - std(temp)/sqrt(length(temp));
        DP2.meanerror{iSubType}(4,iBin) = mean(temp) + std(temp);
        DP2.meanerror{iSubType}(5,iBin) = mean(temp) - std(temp);
    end
end

% Save the plotting results
save([ScratchDir,'ModelType1PlotData'],'DP1','DP2')

disp('Finished Organizing the Data for Plotting')

%% Make the Plot


% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.1;
rmargin = 0.1;

% Set the top and bottom margins in inches
tmargin = 0.23;
bmargin = 0.3;

% Set the horizontal distance between plots in inches
hspace = [0.4,0.8];

% Set the vertical distance between plots in inches
vspace = 0.5;

% Set Paper size
papersize = [7.4 5.25];

% Set the figure panel dimensions (number of horizontal figures by number
% of vertical figures
figdim = [3,4];

% Set the height proportions for each row
heightprop = [1,1,1,1];

% Set the width proportions for each column
widthprop = [1,3,2];

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
AxLabelFS = 6;
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

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');



% Make the diagrams
% Note that in the publication figures these were added in Illustrator

for iST = 1:4
    f = subplot('Position',[LeftCoord(1),BottomCoord(iST),width(1),height(iST)]);
    hold on
    
    if iST == 1
        title('Single Stimulus','FontSize',TitleFS)
    elseif iST == 2
        title('Two Stimuli','FontSize',TitleFS)
    elseif iST == 3
        title('Delayed Stimulus','FontSize',TitleFS)
    elseif iST == 4
        title('Nonlinear Stimulus','FontSize',TitleFS)
    end
    
    axis off
    
end



% Make the example voltage, spike, stimulus plots

for iST = 1:4
    
    f = subplot('Position',[LeftCoord(2),BottomCoord(iST),width(2),height(iST)]);
    hold on
    
    xMin = 0;
    xMax = max(DP1.time{iST});
    yMax = max(DP1.mV{iST});
    yMin = min(DP1.mV{iST});
    dy = yMax - yMin;
    yMax = yMax + 0.1*dy;
    yMin = yMin - 0.1*dy;
    
    % Transform the current to voltage units for plotting
    cMax = max(DP1.current{iST});
    cMin = min(DP1.current{iST});
    dc = cMax - cMin;
    cMax = cMax + 0.05*dc;
    cMin = cMin - 0.05*dc;
    m = (yMax - yMin)/(cMax - cMin);
    b = yMin - m*cMin;
    
    yMax = yMax + 0.2*(yMax - yMin);
    
    % Plot the current, voltage, and spikes
    
    plot(DP1.time{iST},DP1.mV{iST},'k','LineWidth',VoltageLW)
    plot(DP1.time{iST},m*DP1.current{iST} + b,'Color',PColor{2},'LineWidth',VoltageLW)
    for iSpike = 1:length(DP1.spkTimes{iST})
        line([DP1.spkTimes{iST}(iSpike),DP1.spkTimes{iST}(iSpike)],[yMin + 0.85*(yMax - yMin),yMax],'Color',PColor{1},'LineWidth',VoltageLW)
    end
    
    % If this is subtype 3, show the delay
    if iST == 3
        delay = DP1.delay(iST);
        dis = 0.1;
        line([500,500 + delay],[yMin + dis*(yMax - yMin),yMin + dis*(yMax - yMin)],'Color',PColor{4},'LineWidth',DelayLW)
        text(500 + delay + 0.01*(xMax - xMin),yMin + dis*(yMax - yMin),'Delay','FontSize',DelayFS,'Color',PColor{4})
    end
    
    % Add the second y-axis on the right
    line([xMax,xMax],[yMin,yMax],'Color',PColor{2},'LineWidth',AxLW)
    if iST == 4
        cTicks = [0,500,1000];
    else
        cTicks = [0,250,500];
    end
    for iTick = 1:length(cTicks)
        line([xMax - cTickL*(xMax - xMin)/(width(2)*papersize(1)),xMax],[m*cTicks(iTick) + b,m*cTicks(iTick) + b],'Color',PColor{2},'LineWidth',AxLW)
        text(xMax + cTickspacing*(xMax - xMin)/(width(2)*papersize(1)),m*cTicks(iTick) + b,num2str(cTicks(iTick)),'FontSize',UnitFS,'Color',PColor{2})
    end
    if iST == 4
        text(xMax + cTickLabelspacing(2)*(xMax - xMin)/(width(2)*papersize(1)),mean([yMax,yMin]),'Stimulus Current (pA)','Color',PColor{2},...
            'VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-90,'FontSize',AxLabelFS)
    else
        text(xMax + cTickLabelspacing(1)*(xMax - xMin)/(width(2)*papersize(1)),mean([yMax,yMin]),'Stimulus Current (pA)','Color',PColor{2},...
            'VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-90,'FontSize',AxLabelFS)
    end
    
    xlim([xMin,xMax])
    ylim([yMin,yMax])
    set(gca,'FontSize',UnitFS)
    xlabel('Time (ms)','FontSize',AxLabelFS)
    ylabel('Membrane Voltage (mV)','FontSize',AxLabelFS)
    title('Example Voltage/Stimulus Traces','FontSize',TitleFS)
    
end



% Plot the Information Box Plots

for iST = 1:4
    
    f = subplot('Position',[LeftCoord(3),BottomCoord(iST),width(3),height(iST)]);
    hold on
    
    yMax = 1;
    yMin = 0;
    dy = yMax - yMin;
    yMax = yMax + 0.1*dy;
    yMin = yMin - 0.1*dy;
    
    nBins = size(DP2.boxpoints{iST},2);
    
    line([0,0],[yMin,yMax],'Color',PColor{3},'LineStyle','--','LineWidth',StimLineLW)
    text(-15,mean([yMin,yMax]),'Stim On','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',90,'FontSize',StimLabelFS,'Color',PColor{3})
    line([500,500],[yMin,yMax],'Color',PColor{2},'LineStyle','--','LineWidth',StimLineLW)
    text(500 + 15,mean([yMin,yMax]),'Stim Off','VerticalAlignment','bottom','HorizontalAlignment','center','Rotation',-90,'FontSize',StimLabelFS,'Color',PColor{2})
    
    for iBin = 1:nBins
        line([DP2.time{iST}(iBin),DP2.time{iST}(iBin)],[DP2.meanerror{iST}(4,iBin),DP2.meanerror{iST}(5,iBin)],'Color','k','LineWidth',BPThinLW)
        scatter(DP2.time{iST}(iBin),DP2.meanerror{iST}(1,iBin),BPDotSz,'k','filled')
    end
    
    ylim([yMin,yMax])
    xlim([-200,800])
    
    % If this is subtype 3, show the delay
    if iST == 3
        delay = DP2.delay(iST);
        dis = 0.25;
        line([0,delay],[yMin + dis*(yMax - yMin),yMin + dis*(yMax - yMin)],'Color',PColor{4},'LineWidth',DelayLW)
        text(delay + 10,yMin + dis*(yMax - yMin),'Delay','FontSize',DelayFS,'Color',PColor{4})
    end
    
    set(gca,'FontSize',UnitFS)
    xlabel('Time (ms)','FontSize',AxLabelFS)
    ylabel('Information (bits)','FontSize',AxLabelFS)
    title('Spike Count/Stimulus Encoding','FontSize',TitleFS)
    
end


% Print the Figure
cd(ScratchDir(1:(end - 1)))
print(F1,'-dpdf','-painters','neuroDemo1Panel')
cd(SimDir)