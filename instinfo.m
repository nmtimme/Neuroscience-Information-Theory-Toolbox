%% INSTINFO - instantaneous information analysis
% Calculates various information theoretic quantities for trial based or
% single trial data. Single trial data are analyzed through time, while
% trial based data are analyzed across trials at each time bin.
%
% Syntax: [InfoVal,p,s,MCInfoVals] = instinfo(StatesRaster, Method, VariableIDs, varargin)
%
% Input:
%   StatesRaster (cell array or double array) - trial state data. If a
%     double array is used, it should be number of variables by number of
%     time bins by number of trials. Each element should be an integer
%     state number (state number = 1, 2, 3, ...). If a cell array is used,
%     it should have only one dimension and each element should be a double
%     array with the dimensions listed above. Each element of the cell
%     array is referred to as a 'data category'.
%   Method (string) - sets the type of information measure to apply. (See
%     option sets below)
%   VariableIDs (cell array) - sets the variables used in the information
%     measure. The array is generally N by 3 where N is the number of 
%     variables. The first column contains an integer that sets the data 
%     category. If there is only one data category, this column can be 
%     removed. The second column contains the variable number in the data 
%     category. This number corresponds to the first subscript in the 
%     StatesRaster array. The third column contains the time bin value. 
%     This corresponds to the second subscript in the StatesRaster array. 
%     (See option sets below.)
%     *Side Note: For single trial data, the precise values of the time
%       bins in the third column are not relevant, only the relative values
%       will affect the analysis. The relative values will set delays
%       between the states used. For instance, for Pairwise Mutual
%       Information, {1,5,2;1,1,3} is equivalent to {1,5,5;1,1,6} because 
%       3 - 2 = 6 - 5, though different from {1,5,2;1,1,4} because 
%       3 - 2 ~= 4 - 2. Positive integers must be used. Frequently, it is
%       best to set the earliest variable as time 1 and proceed from there.
%
% Input Option Sets
%   Entropy - measures the entropy of a given source. The source to be used
%     is listed in the VariableIDs. The Monte Carlo options do not affect 
%     entropy calculations. InfoVal output will be a single number in units
%     of bits.
%     Example 1
%       Method = 'Ent'
%       VariableIDs = {1,7,3}
%       Explanation: These inputs will produce an entropy calculation for
%       data category 1 variable 7 at time bin 3.
%     Example 2
%       Method = 'Ent'
%       VariableIDs = {3,2,5}
%       Explanation: These inputs will produce an entropy calculation for
%       data category 3 variable 2 at time bin 5.
%   Joint Entropy - measures the entropy of several sources jointly. The
%     sources to be used are listed in the VariableIDs. The Monte Carlo
%     options do not affect entropy calculation. InfoVal output will be a
%     single number in units of bits.
%     Example 1
%       Method = 'JointEnt'
%       VariableIDs = {1,7,3;...
%                      2,5,6}
%       Explanation: These inputs will produce a joint entropy calculation 
%       for data category 1 variable 7 at time bin 3 and data category 2
%       variable 5 at time bin 6.
%     Example 2
%       Method = 'JointEnt'
%       VariableIDs = {1,7,3;...
%                      2,5,6;...
%                      2,9,10}
%       Explanation: These inputs will produce a joint entropy calculation 
%       for data category 1 variable 7 at time bin 3, data category 2
%       variable 5 at time bin 6, and data category 2 variable 9 at time
%       bin 10.
%   Joint Conditional Entropy - measures the conditional entropy between any
%     two groups of sources, each taken jointly. The sources are listed in 
%     VariableIDs. The grouping is indicated using the method name 
%     ('CondEnti' where the first i variables listed in VariableIDs will be
%     considered as one joint set and conditioned on the remaining variable
%     considered as another joint set, with i being less than or equal to
%     9). The Monte Carlo options do not affect entropy calculation. 
%     InfoVal output will be a single number in units of bits.
%     Example 1
%       Method = 'CondEnt1'
%       VariableIDs = {1,7,3;...
%                      2,5,6}
%       Explanation: These inputs will produce a conditional entropy 
%       calculation with data category 1 variable 7 at time bin 3 
%       conditioned on data category 2 variable 5 at time bin 6.
%     Example 2
%       Method = 'CondEnt2'
%       VariableIDs = {1,7,3;...
%                      2,5,6;...
%                      2,9,10}
%       Explanation: These inputs will produce a conditonal entropy 
%       calculation with data category 1 variable 7 at time bin 3 jointly
%       with data category 2 variable 5 at time bin 6 on data category 2
%       variable 9 at time bin 10. 
%   Pairwise Mutual Information - measures the pairwise mutual information 
%     between any two sources. The sources are listed in VariableIDs. If 
%     the Monte Carlo option is selected, the first variable will be 
%     randomized, though this is equivalent to randomizing the second 
%     variable because the joint distribution is bivariate. InfoVal output 
%     will be a single number in units of bits.
%     Example 1
%       Method = 'PairMI'
%       VariableIDs = {1,6,5;...
%                      1,10,5}
%       Explanation: These inputs will produce a pairwise mutual
%       information measurement between data category 1 variables 6 and 10 
%       contemporaneously at time bin 5.
%     Example 2
%       Method = 'PairMI'
%       VariableIDs = {8,5,2;...
%                      9,1,3}
%       Explanation: These inputs will produce a pairwise mutual
%       information measurement between data category 8 variable 5 at time 
%       bin 2 and data category 9 variable 1 at time bin 3. 
%   Joint Mutual Information - measures the mutual information between any
%     two groups of sources, each taken jointly. The sources are listed in 
%     VariableIDs. The grouping is indicated using the method name 
%     ('JointMIi' where the first i variables listed in VariableIDs will be
%     considered as one joint set and the remaining variables will be 
%     considered as the other joint set, with i being less than or equal to
%     9). If the Monte Carlo option is selected, the first joint set of 
%     variables will be randomized together, though this is equivalent to 
%     randomizing the second joint set of variables because the 
%     distribution is bivariate. InfoVal output will be a single number in 
%     units of bits. 
%     Example 1
%       Method = 'JointMI2'
%       VariableIDs = {1,6,5;...
%                      4,10,5;...
%                      9,2,3}
%       Explanation: These inputs will produce a mutual information
%       measurement between the joint set of data category 1 variable 6 at
%       time bin 5 with data category 4 variable 10 at time bin 5, and data
%       category 9 variable 2 at time bin 3.
%     Example 2
%       Method = 'JointMI1'
%       VariableIDs = {6,5,2;...
%                      7,1,3;...
%                      8,4,10}
%       Explanation: These inputs will produce a mutual information
%       measurement between data category 6 variable 5 at time bin 2, and
%       the joint set of data category 7 variable 1 at time bin 3 with data
%       category 8 variable 4 at time bin 10. 
%   Transfer Entropy - measures the transfer entropy between any two
%     sources. The first variable listed in the VariableIDs will be the 
%     receiving variable in the future, the second variable listed will be 
%     the receiving variable in the past, and the third variable listed 
%     will be the transmitting variable in the past. Typically, the time 
%     bin for both past variables are prior to the future state (often one 
%     time bin earlier). An error will be issued if the times for the past 
%     states are later than or equal to the time for the future state. If 
%     the Monte Carlo option is selected, only the transmitter variable 
%     will be randomized to preserve the autoprediction in the receiver. 
%     InfoVal output will be a single number in units of bits.
%     Example 1
%       Method = 'TE'
%       VariableIDs = {1,3,7;...
%                      1,3,6;...
%                      2,8,6}
%       Explanation: These inputs will calculate TE with the transmitter
%       variable being data category 2 variable 8 at time bin 6, the 
%       receiver past state will be data category 1 variable 3 at time bin 
%       6, and the receiver future state will be data category 1 variable 3
%       at time bin 7.
%     Example 2
%       Method = 'TE'
%       VariableIDs = {3,5,7;...
%                      3,5,6;...
%                      1,2,4}
%       Explanation: These inputs will calculate TE with the transmitter
%       variable being data category 1 variable 2 at time bin 4, the 
%       receiver past state will be data category 3 variable 5 at time bin 
%       6, and the receiver future state will be data category 3 variable 5
%       at time bin 7.
%   2-Variable PID - measures the 2 X variable, 1 Y variable partial
%     information terms. The first variable in VariableIDs will be the Y 
%     variable. The second and third variables in VariableIDs will be the 
%     first and second X variables, respectively. If the Monte Carlo option
%     is selected, both of the X variables will be randomized. InfoVal 
%     output will contain four values (redundancy, unique X1, unique X2, 
%     and synergy). The p output will contain 4 p-values (corresponding to 
%     the terms in InfoVal). 
%     Example 1
%       Method = '2PID'
%       VariableIDs = {1,12,3;...
%                      4,8,3;...
%                      4,2,3}
%       Explanation: These inputs will calculate the 2-variable PID with
%       data category 1 variable 12 at time bin 3 as the Y variable, data 
%       category 4 variable 8 at time bin 3 as the X1 variable, and data 
%       category 4 variable 2 at time bin 3 as the X2 variable.
%     Example 2
%       Method = '2PID'
%       VariableIDs = {2,1,5;...
%                      8,3,7;...
%                      9,4,6}
%       Explanation: These inputs will calculate the 2-variable PID with
%       data category 2 variable 1 at time bin 5 as the Y variable, data 
%       category 8 variable 3 at time bin 7 as the X1 variable, and data 
%       category 9 variable 4 at time bin 6 as the X2 variable.
%   3-Variable PID - measures the 3 X variable, 1 Y variable partial
%     information terms. The first variable in VariableIDs will be the Y 
%     variable. The second,  third, and fourth variables in VariableIDs 
%     will be the first, second, and third X variables, respectively. If 
%     the Monte Carlo option is  selected, all three X variables will be 
%     randomized. InfoVal output will contain 18 values (see PIDTermLabels 
%     for term identity). The p output will contain 18 p-values 
%     (corresponding to the terms in InfoVal).
%     Example 1
%       Method = '3PID'
%       VariableIDs = {1,12,5;...
%                      3,8,5;...
%                      3,2,5;...
%                      4,7,5}
%       Explanation: These inputs will calculate the 3-variable PID with
%       data category 1 variable 12 at time bin 5 as the Y variable, data 
%       category 3 variable 8 at time bin 5 as the X1 variable, data 
%       category 3 variable 2 at time bin 5 as the X2 variable, and data 
%       category 4 variable 7 at time bin 5 as the X3 variable.
%     Example 2
%       Method = '3PID'
%       VariableIDs = {9,1,5;...
%                      10,2,6;...
%                      11,3,7;...
%                      12,4,8}
%       Explanation: These inputs will calculate the 2-variable PID with
%       data category 9 variable variable 1 at time bin 5 as the Y 
%       variable, data category 10 variable 2 at time bin 6 as the X1 
%       variable, data category 11 variable 3 at time bin 7 as the X2 
%       variable, and data category 12 variable 4 at time bin 8 as the X3 
%       variable.
%   3-Variable TE - measures the multivariate transfer entropy from 2
%     sources to a third source. The first variable listed in the
%     VariableIDs will be the receiving variable in the future. The second
%     variable listed in the VariableIDs will be the receiving variable in
%     the past. The third variable listed in the VariablesIDs will the
%     first transmitting variable in the past. The fourth variable listed
%     in the VariablesIDs will be the second transmitting variable in the
%     past. If any of the past variables have times later than or equal to
%     the future state, an error will be issued. If the Monte Carlo option 
%     is selected, only the transmitter variables will be randomized to
%     preserve the autoprediction in the receiver. InfoVal output will have
%     4 values (redundancy, unique transmitter 1, unique transmitter 2,
%     synergy), each in units of bits. The p output will have 4 p-values
%     corresponding to each InfoVal.
%     Example 1
%       Method = '3TE'
%       VariableIDs = {2,1,10;...
%                      2,1,9;...
%                      3,5,9;...
%                      4,6,9}
%       Explanation: These inputs will calculate the 3-variable TE with
%       data category 2 variable 1 at time bin 10 as the future state of
%       the receiver, data category 2 variable 1 at time bin 9 as the past 
%       state of the receiver, data category 3 variable 5 at time bin 9 as 
%       the state of transmitter 1, and data category 4 variable 6 at time 
%       bin 9 as the state of transmitter 2.
%     Example 2
%       Method = '3TE'
%       VariableIDs = {3,1,10;...
%                      3,1,9;...
%                      4,5,7;...
%                      4,2,7}
%       Explanation: These inputs will calculate the 3-variable TE with
%       data category 3 variable 1 at time bin 10 as the future state of
%       the receiver, data category 3 variable 1 at time bin 9 as the past 
%       state of the receiver, data category 4 variable 5 at time bin 7 as 
%       the state of transmitter 1, and data category 4 variable 2 at time 
%       bin 7 as the state of transmitter 2.
%   Information Gain - measures the information gained by one variable
%     (call it the receiver) about another variable (call it the signal) 
%     beyond the information present in the receiver about the signal in 
%     the past. The first variable listed in VariableIDs will be the 
%     receiver in the future. The second variable listed in VariableIDs
%     will be the receiver in the past. The third variable listed in
%     VariableIDs will be the signal. An error will be produced if the past
%     state is after or at the same time as the future state. If the Monte 
%     Carlo option is used, the signal states are randomized. InfoVal will 
%     have one value measured in bits.
%     Example 1
%       Method = 'InfoGain'
%       VariableIDs = {2,4,6;...
%                      2,4,5;...
%                      3,1,5}
%       Explanation: These inputs will calculate the information gained by
%       data category 2 variable 4 at time bin 6 about data category 3 
%       variable 1 at time bin 5 beyond the information present in data 
%       category 2 variable 4 at time bin 5 about data category 3 variable 
%       1 at time bin 5.
%     Example 2
%       Method = 'InfoGain'
%       VariableIDs = {7,4,6;...
%                      7,4,5;...
%                      2,1,3}
%       Explanation: These inputs will calculate the information gained by
%       data category 7 variable 4 at time bin 6 about data category 2 
%       variable 1 at time bin 3 beyond the information present in data 
%       category 7 variable 4 at time bin 5 about data category 2 variable 
%       1 at time bin 3.
%   Information Transmission - measures the information transmitted from
%     one source (call it the transmitter) to another source (call it the
%     receiver) about a third source (call it the signal). The first
%     variable listed in VariableIDs is the receiver in the future, the
%     second is the receiver in the past, the third is the signal, and the
%     fourth is the transmitter in the past. An error will be issued if
%     the past states are after the future states. If the Monte Carlo 
%     option is used, the signal states are randomized. InfoVal will have 
%     one value measured in bits.
%     Example 1
%       Method = 'InfoTrans'
%       VariableIDs = {2,4,6;...
%                      2,4,5;...
%                      3,1,5;...
%                      8,12,5}
%       Explanation: These inputs will calculate the information
%       transmitted from data category 8 variable 12 at time bin 5 to data 
%       category 2 variable 4 at time bin 6 about data category 3 
%       variable 1 at time bin 5.
%     Example 2
%       Method = 'InfoTrans'
%       VariableIDs = {2,4,6;...
%                      2,4,5;...
%                      9,1,3;...
%                      8,12,3}
%       Explanation: These inputs will calculate the information
%       transmitted from data category 8 variable 12 at time bin 3 to data 
%       category 2 variable 4 at time bin 6 about data category 9 variable 
%       1 at time bin 3.
%   
%     
%
% Variable Inputs:
%   (..., 'MCOpt', MCOpt) - specifies whether to use the Monte Carlo method
%     to measure a p-value (MCOpt = 'on') or not (MCOpt = 'off') (string)
%     (default: 'off')
%   (..., 'MCpThresh', MCpThresh) - sets the p-value cutoff. If the
%     algorithm finds that the p-value will be above the threshold, the
%     calculation ceases (scalar double) (default: 0.001)
%   (..., 'MCnSamples', MCnSamples) - sets the number of Monte Carlo trials
%     to run (scalar double) (default: 5000)
%   (..., 'nullModel', nullModel) - allows the user to use a previously
%     established null model probability mass function to calculate
%     significance instead of Monte Carlo. This can be especially efficient
%     when many information theory analyses have the same underlying null
%     model. Note that this option is only currently available for 'PairMI'
%     and '2PID'. (structure with fields:
%       nullModel.Method - the information theory analysis method (see
%         above) used to generate the null distribution (string). This must
%         match Method, otherwise the function will employ Monte Carlo, if
%         requested.
%       nullModel.rFact - a rounding factor used to correct for rounding
%         errors near the system resolution (double). The information value
%         from the real data will be rounded to the scale set by rFact for
%         comparison to the null model pmf.
%       nullModel.res - the resolution of the null model (double). This
%         value sets the lowest possible p-value for the null model.
%       nullModel.nCounts - the number of observations for each state of
%         each variable listed in VariableIDs (in order) (cell array with
%         each cell as a column vector with state count numbers). This
%         information will be compared against the data input into instinfo
%         to ensure that the null models are identical. If the null model
%         cannot be used, Monte Carlo will be employed, if requested.
%       nullModel.info - the unique information values produced by the null
%         model (cell array with each cell as a vector of ascending
%         information values order matched to the information value outputs
%         of instinfo). These values serve as the independent variable in
%         the probability mass function. If only one information value is
%         supplied, nullModel.info can be a vector instead of a cell.
%       nullModel.pmf - the probability of each unique information value
%         from the corresponding location in nullModel.info. These values
%         serve as the dependent variable in the probability mass function.
%         The p-value will be 1 - the cumulative mass function at the
%         information value observed in the real data. The type (cell vs.
%         vector) and size of nullModel.pmf must match nullModel.info.
%
% Outputs:
%   InfoVal (scalar or vector double) - information value in bits (see
%     option set examples above for precise meaning). If multiple
%     information values are output (e.g., Method = '2PID'), then InfoVal
%     will be a column vector.
%   p (scalar or vector double) - p-value from the Monte Carlo estimation
%     (will be NaN if the Monte Carlo option is not selected). If multiple
%     information values are output (e.g., Method = '2PID'), then p will be
%     a column vector corresponding to InfoVal. Note that the minimum
%     p-value is set by the resolution from the supplied null model or the
%     number of Monte Carlo trials (1/MCnSamples). In both cases, p-value
%     results of 0 are reset to half the resolution. 
%   MCInfoVals (vector or array double) - information values from the Monte
%     Carlo trials (will be NaN if the Monte Carlo option is not selected).
%     If multiple information values are output (e.g., Method = '2PID'),
%     then MCInfoVals will be a number of information values by number of
%     Monte Carlo trials matrix.
%   s (double) - a value that indicates whether the Monte Carlo approach
%     was used to assess significance (s = 0), the supplied null model was
%     used (s = 1), or no significance testing was applied (s = NaN).
%
%
% Other m-files required: EntropyY, MutualInfo, TE2, PID
% Subfunctions: none
% MAT-files required: 2PIDMats.mat, 3PIDMats.mat, TE3Redux.mat
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% August 2014; Last revision: 12-May-2016


function [InfoVal,p,s,MCInfoVals] = instinfo(StatesRaster, Method, VariableIDs, varargin)
%% Parse command line for parameters
MCOpt = 'off';
MCpThresh = 0.001;
MCnSamples = 5000;
p = NaN;
MCInfoVals = NaN;
NullModelOpt = 'off';
s = NaN;

iVarArg = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'MCOpt',        MCOpt = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'MCpThresh',    MCpThresh = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'MCnSamples',   MCnSamples = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'nullModel',    NullModelOpt = 'on'; nullModel = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(INSTINFO) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

%% Perform Initial Calculations

% Error check inputs and possible reorganize to ease later processing
if iscell(StatesRaster)
    if length(StatesRaster) ~= 1
        if size(VariableIDs,2) ~= 3
            error('VariableIDs is not the appropriate size. Probably the data category is missing.')
        end
    else
        if size(VariableIDs,2) == 2
            NewVariableIDs = cell([size(VariableIDs,1),3]);
            for iVar = 1:size(VariableIDs,1)
                NewVariableIDs{iVar,1} = 1;
                NewVariableIDs{iVar,2} = VariableIDs{iVar,1};
                NewVariableIDs{iVar,3} = VariableIDs{iVar,2};
            end
            VariableIDs = NewVariableIDs;
        elseif size(VariableIDs,2) ~= 3
            error('VariableIDs is not the appropriate size.')
        end
    end
else
    StatesRaster = {StatesRaster};
    if size(VariableIDs,2) == 2
        NewVariableIDs = cell([size(VariableIDs,1),3]);
        for iVar = 1:size(VariableIDs,1)
            NewVariableIDs{iVar,1} = 1;
            NewVariableIDs{iVar,2} = VariableIDs{iVar,1};
            NewVariableIDs{iVar,3} = VariableIDs{iVar,2};
        end
        VariableIDs = NewVariableIDs;
    end
end
    
% Detect single trial data
STTest = zeros([1,size(VariableIDs,1)]);
for iVar = 1:size(VariableIDs,1)
    if size(StatesRaster{VariableIDs{iVar,1}},3) == 1
        STTest(iVar) = 1;
    end
end
if sum(STTest) == size(VariableIDs,1)
    STFlag = true;
    STLength = zeros([1,size(VariableIDs,1)]);
    for iVar = 1:size(VariableIDs,1)
        STLength(iVar) = size(StatesRaster{VariableIDs{iVar,1}},2);
    end
    if ~all(STLength == STLength(1))
        error('All recording lengths must be identical for single trial data.')
    end
elseif sum(STTest) == 0
    STFlag = false;
    nTrials = zeros([1,size(VariableIDs,1)]);
    for iVar = 1:size(VariableIDs,1)
        nTrials(iVar) = size(StatesRaster{VariableIDs{iVar,1}},3);
    end
    if ~all(nTrials == nTrials(1))
        error('All data must have the same number of trials for trial based data.')
    end
else
    error('Mixture of single trial and trial based data is not allowed.')
end
        

%% Perform the Relevant Information Measurement

if strcmp(Method,'Ent')
    
    %% Perform the Entropy Calculation
    
    % Parse the inputs into trial states
    if STFlag
        x = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},:));
    else
        x = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},VariableIDs{1,3},:));
    end
    
    % Error check the integer states
    if ~isequal(unique(x),(1:length(unique(x))))
        [B,I,x] = unique(x);
    end
    
    % Make a counts matrix
    Counts = accumarray({x},ones(size(x)),[length(unique(x)),1]);
    
    % Calculate the entropy
    InfoVal = EntropyY(Counts);
    
elseif strcmp(Method,'JointEnt')
    
    %% Perform the Joint Entropy Calculation
    
    nJoint = length(VariableIDs(:,3));
    if STFlag
        t1 = cell2mat(VariableIDs(:,3)) - min(cell2mat(VariableIDs(:,3)));
        t2 = cell2mat(VariableIDs(:,3)) - max(cell2mat(VariableIDs(:,3)));
        nT = size(StatesRaster{VariableIDs{1,1}},2) - max(t1);
        x = NaN([nT,nJoint]);
        
        % Parse the inputs into trial states
        for iVar = 1:nJoint
            x(:,iVar) = squeeze(StatesRaster{VariableIDs{iVar,1}}(VariableIDs{iVar,2},(t1(iVar) + 1):(end + t2(iVar))));
        end
    else
        nTrials = size(StatesRaster{VariableIDs{1,1}},3);
        x = NaN([nTrials,nJoint]);
        
        % Parse the inputs into trial states
        for iVar = 1:nJoint
            x(:,iVar) = squeeze(StatesRaster{VariableIDs{iVar,1}}(VariableIDs{iVar,2},VariableIDs{iVar,3},:));
        end
    end
    [B,I,x] = unique(x,'rows');
    
    % Error check the integer states
    if ~isequal(unique(x),(1:length(unique(x))))
        [B,I,x] = unique(x);
    end
    
    % Make a counts matrix
    Counts = accumarray({x},ones(size(x)),[length(unique(x)),1]);
    
    % Calculate the entropy
    InfoVal = EntropyY(Counts);
    
elseif strcmp(Method(1:(end - 1)),'CondEnt')
    
    %% Perform the Conditional Entropy Calculation
    
    nJoint1 = str2double(Method(end));
    nJoint2 = size(VariableIDs,1) - nJoint1;
    if nJoint2 < 1
        error('VariableIDs number of neurons does not match the method label.')
    end
    if STFlag
        t1 = cell2mat(VariableIDs(:,3)) - min(cell2mat(VariableIDs(:,3)));
        t2 = cell2mat(VariableIDs(:,3)) - max(cell2mat(VariableIDs(:,3)));
        nT = size(StatesRaster{VariableIDs{1,1}},2) - max(t1);
        x = NaN([nT,nJoint1]);
        y = NaN([nT,nJoint2]);
        
        % Parse the inputs into trial states
        for iVar = 1:nJoint1
            x(:,iVar) = squeeze(StatesRaster{VariableIDs{iVar,1}}(VariableIDs{iVar,2},(t1(iVar) + 1):(end + t2(iVar))));
        end
        for iVar = (nJoint1 + 1):(nJoint1 + nJoint2)
            y(:,iVar - nJoint1) = squeeze(StatesRaster{VariableIDs{iVar,1}}(VariableIDs{iVar,2},(t1(iVar) + 1):(end + t2(iVar))));
        end
    else
        nTrials = size(StatesRaster{VariableIDs{1,1}},3);
        x = NaN([nTrials,nJoint1]);
        y = NaN([nTrials,nJoint2]);
        
        % Parse the inputs into trial states
        for iVar = 1:nJoint1
            x(:,iVar) = squeeze(StatesRaster{VariableIDs{iVar,1}}(VariableIDs{iVar,2},VariableIDs{iVar,3},:));
        end
        for iVar = (nJoint1 + 1):(nJoint1 + nJoint2)
            y(:,iVar - nJoint1) = squeeze(StatesRaster{VariableIDs{iVar,1}}(VariableIDs{iVar,2},VariableIDs{iVar,3},:));
        end
    end
    [B,I,x] = unique(x,'rows');
    [B,I,y] = unique(y,'rows');
    
    % Error check the integer states
    if ~isequal(unique(x),(1:length(unique(x))))
        [B,I,x] = unique(x);
    end
    if ~isequal(unique(y),(1:length(unique(y))))
        [B,I,y] = unique(y);
    end
    
    % Make the counts matrices
    CountsXY = accumarray({x,y},ones(size(x)),[length(unique(x)),length(unique(y))]);
    CountsX = accumarray({x},ones(size(x)),[length(unique(x)),1]);
    
    % Calculate the Conditional Entropy
    InfoVal = EntropyY(CountsX) - MutualInfo(CountsXY);
    
    % Correct for rounding errors
    if abs(InfoVal) < (10*eps)
        InfoVal = 0;
    end

elseif strcmp(Method,'PairMI')
    
    %% Perform the Pairwise Mutual Information Calculation
    
    % Parse the inputs into trial states
    if STFlag
        t1 = cell2mat(VariableIDs(:,3)) - min(cell2mat(VariableIDs(:,3)));
        t2 = cell2mat(VariableIDs(:,3)) - max(cell2mat(VariableIDs(:,3)));
        x = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},(t1(1) + 1):(end + t2(1))));
        y = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},(t1(2) + 1):(end + t2(2))));
    else
        x = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},VariableIDs{1,3},:));
        y = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},VariableIDs{2,3},:));
    end
    
    % Error check the integer states
    if ~isequal(unique(x),(1:length(unique(x))))
        [B,I,x] = unique(x);
    end
    if ~isequal(unique(y),(1:length(unique(y))))
        [B,I,y] = unique(y);
    end
    
    % Make a counts matrix
    Counts = accumarray({x,y},ones(size(x)),[length(unique(x)),length(unique(y))]);
    
    % Calculate the pairwise mutual information
    InfoVal = MutualInfo(Counts);
    
    % If the user supplied a null model pmf, use it
    if strcmp(NullModelOpt,'on')
        try
            % Try to use the null model, but allow for the possibility that
            % the null model does not apply to these data, which will
            % produce an error.
            p = nullpmf(InfoVal,Counts,nullModel,Method);
            
            % Record that the null model was used
            s = 1;
        catch
            % If the null model does not apply to these data, turn the
            % option off so the Monte Carlo can run, if necessary. Also,
            % warn the user
            NullModelOpt = 'off';
            warning('The supplied null model does not apply to these data.')
        end
    end
    
    % If the user requested Monte Carlo, perform it
    if strcmp(MCOpt,'on') && strcmp(NullModelOpt,'off')
        
        % Figure out how many failures we're allowed
        nFailsThresh = ceil(MCpThresh * MCnSamples);
        iSample = 1;
        iFails = 0;
        MCInfoVals = NaN([1,MCnSamples]);
        nx = length(unique(x));
        ny = length(unique(y));
        lx = length(x);
        
        % Perform the Monte Carlo Trials
        while (iFails < nFailsThresh) && (iSample <= MCnSamples)
            % Counts = accumarray({x(randperm(lx)),y},ones(size(x)),[nx,ny]);
            MCInfoVals(iSample) = MutualInfo(accumarray({x(randperm(lx)),y},ones(size(x)),[nx,ny]));
            iFails(MCInfoVals(iSample) >= InfoVal) = iFails(MCInfoVals(iSample) >= InfoVal) + 1;
            iSample = iSample + 1;
        end
        
        % Calculate the p-value
        p = iFails / (iSample - 1);
        
        % Correct for the resolution of the Monte Carlo trials
        p(p == 0) = 1/(2*MCnSamples);
        
        % Record that the Monte Carlo method was used
        s = 0;
    end
    
elseif strcmp(Method(1:(end - 1)),'JointMI')
    
    %% Perform the Joint Mutual Information Calculation
    
    nJoint1 = str2double(Method(end));
    nJoint2 = size(VariableIDs,1) - nJoint1;
    if nJoint2 < 1
        error('VariableIDs number of neurons does not match the method label.')
    end
    if STFlag
        t1 = cell2mat(VariableIDs(:,3)) - min(cell2mat(VariableIDs(:,3)));
        t2 = cell2mat(VariableIDs(:,3)) - max(cell2mat(VariableIDs(:,3)));
        nT = size(StatesRaster{VariableIDs{1,1}},2) - max(t1);
        x = NaN([nT,nJoint1]);
        y = NaN([nT,nJoint2]);
        
        % Parse the inputs into trial states
        for iVar = 1:nJoint1
            x(:,iVar) = squeeze(StatesRaster{VariableIDs{iVar,1}}(VariableIDs{iVar,2},(t1(iVar) + 1):(end + t2(iVar))));
        end
        for iVar = (nJoint1 + 1):(nJoint1 + nJoint2)
            y(:,iVar - nJoint1) = squeeze(StatesRaster{VariableIDs{iVar,1}}(VariableIDs{iVar,2},(t1(iVar) + 1):(end + t2(iVar))));
        end
    else
        nTrials = size(StatesRaster{VariableIDs{1,1}},3);
        x = NaN([nTrials,nJoint1]);
        y = NaN([nTrials,nJoint2]);
        
        % Parse the inputs into trial states
        for iVar = 1:nJoint1
            x(:,iVar) = squeeze(StatesRaster{VariableIDs{iVar,1}}(VariableIDs{iVar,2},VariableIDs{iVar,3},:));
        end
        for iVar = (nJoint1 + 1):(nJoint1 + nJoint2)
            y(:,iVar - nJoint1) = squeeze(StatesRaster{VariableIDs{iVar,1}}(VariableIDs{iVar,2},VariableIDs{iVar,3},:));
        end
    end
    [B,I,x] = unique(x,'rows');
    [B,I,y] = unique(y,'rows');
    
    % Error check the integer states
    if ~isequal(unique(x),(1:length(unique(x))))
        [B,I,x] = unique(x);
    end
    if ~isequal(unique(y),(1:length(unique(y))))
        [B,I,y] = unique(y);
    end
    
    % Make a counts matrix
    Counts = accumarray({x,y},ones(size(x)),[length(unique(x)),length(unique(y))]);
    
    % Calculate the pairwise mutual information
    InfoVal = MutualInfo(Counts);
    
    % If the user requested Monte Carlo, perform it
    if strcmp(MCOpt,'on')
        
        % Figure out how many failures we're allowed
        nFailsThresh = ceil(MCpThresh * MCnSamples);
        iSample = 1;
        iFails = 0;
        MCInfoVals = NaN([1,MCnSamples]);
        nx = length(unique(x));
        ny = length(unique(y));
        lx = length(x);
        
        % Perform the Monte Carlo Trials
        while (iFails < nFailsThresh) && (iSample <= MCnSamples)
            % Counts = accumarray({x(randperm(lx)),y},ones(size(x)),[nx,ny]);
            MCInfoVals(iSample) = MutualInfo(accumarray({x(randperm(lx)),y},ones(size(x)),[nx,ny]));
            iFails(MCInfoVals(iSample) >= InfoVal) = iFails(MCInfoVals(iSample) >= InfoVal) + 1;
            iSample = iSample + 1;
        end
        
        % Calculate the p-value
        p = iFails / (iSample - 1);
        
        % Correct for the resolution of the Monte Carlo trials
        p(p == 0) = 1/(2*MCnSamples);
        
        % Record that the Monte Carlo method was used
        s = 0;
    end
    
elseif strcmp(Method,'TE')
    
    %% Perform the Transfer Entropy Calculation
    
    % Error check the times and the variable identities
    if VariableIDs{2,3} >= VariableIDs{1,3}
        error('Past state of receiver after the future state of the receiver!')
    elseif VariableIDs{3,3} >= VariableIDs{1,3}
        error('Past state of transmitter after the future state of the receiver!')
    elseif ~isequal(cell2mat(VariableIDs(1,1:2)),cell2mat(VariableIDs(2,1:2)))
        error('Past and future states of the receiver come from different variables!')
    end
    
    % Parse the inputs into trial states
    if STFlag
        t1 = cell2mat(VariableIDs(:,3)) - min(cell2mat(VariableIDs(:,3)));
        t2 = cell2mat(VariableIDs(:,3)) - max(cell2mat(VariableIDs(:,3)));
        x = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},(t1(3) + 1):(end + t2(3))));
        yPast = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},(t1(2) + 1):(end + t2(2))));
        yFuture = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},(t1(1) + 1):(end + t2(1))));
    else
        x = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},VariableIDs{3,3},:));
        yPast = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},VariableIDs{2,3},:));
        yFuture = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},VariableIDs{1,3},:));
    end
    
    % Error check the integer states
    if ~isequal(unique(x),(1:length(unique(x))))
        [B,I,x] = unique(x);
    end
    if ~isequal(unique(yPast),(1:length(unique(yPast))))
        [B,I,yPast] = unique(yPast);
    end
    if ~isequal(unique(yFuture),(1:length(unique(yFuture))))
        [B,I,yFuture] = unique(yFuture);
    end
    
    % Make a counts matrix
    Counts = accumarray({yFuture,yPast,x},ones(size(x)),[length(unique(yFuture)),length(unique(yPast)),length(unique(x))]);
    
    % Calculate the pairwise mutual information
    InfoVal = TE2(Counts);
    
    % If the user requested Monte Carlo, perform it
    if strcmp(MCOpt,'on')
        
        % Figure out how many failures we're allowed
        nFailsThresh = ceil(MCpThresh * MCnSamples);
        iSample = 1;
        iFails = 0;
        MCInfoVals = NaN([1,MCnSamples]);
        nyFuture = length(unique(yFuture));
        nyPast = length(unique(yPast));
        nx = length(unique(x));
        lx = length(x);
        
        % Perform the Monte Carlo Trials
        while (iFails < nFailsThresh) && (iSample <= MCnSamples)
            % Counts = accumarray({yFuture,yPast,x(randperm(lx))},ones(size(x)),[nyFuture,nyPast,nx]);
            MCInfoVals(iSample) = TE2(accumarray({yFuture,yPast,x(randperm(lx))},ones([lx,1]),[nyFuture,nyPast,nx]));
            iFails(MCInfoVals(iSample) >= InfoVal) = iFails(MCInfoVals(iSample) >= InfoVal) + 1;
            iSample = iSample + 1;
        end
        
        % Calculate the p-value
        p = iFails / (iSample - 1);
        
        % Correct for the resolution of the Monte Carlo trials
        p(p == 0) = 1/(2*MCnSamples);
        
        % Record that the Monte Carlo method was used
        s = 0;
    end
    
    
elseif strcmp(Method,'2PID')
    
    %% Perform the 2-Variable PID Calculation
    
    % Parse the inputs into trial states
    if STFlag
        t1 = cell2mat(VariableIDs(:,3)) - min(cell2mat(VariableIDs(:,3)));
        t2 = cell2mat(VariableIDs(:,3)) - max(cell2mat(VariableIDs(:,3)));
        y = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},(t1(1) + 1):(end + t2(1))));
        x1 = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},(t1(2) + 1):(end + t2(2))));
        x2 = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},(t1(3) + 1):(end + t2(3))));
    else
        y = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},VariableIDs{1,3},:));
        x1 = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},VariableIDs{2,3},:));
        x2 = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},VariableIDs{3,3},:));
    end
    
    % Error check the integer states
    if ~isequal(unique(x1),(1:length(unique(x1))))
        [B,I,x1] = unique(x1);
    end
    if ~isequal(unique(x2),(1:length(unique(x2))))
        [B,I,x2] = unique(x2);
    end
    if ~isequal(unique(y),(1:length(unique(y))))
        [B,I,y] = unique(y);
    end
    
    % Make a counts matrix
    Counts = accumarray({y,x1,x2},ones(size(y)),[length(unique(y)),length(unique(x1)),length(unique(x2))]);
    
    % Load the necessary matrices
    load('2PIDMats.mat')
    
    % Calculate the PID values
    InfoVal = PID(Counts,SourcesMat,SetsMat,TransMat);
    
    % If the user supplied a null model pmf, use it
    if strcmp(NullModelOpt,'on')
        try
            % Try to use the null model, but allow for the possibility that
            % the null model does not apply to these data, which will
            % produce an error.
            p = nullpmf(InfoVal,Counts,nullModel,Method);
            
            % Record that the null model was used
            s = 1;
        catch
            % If the null model does not apply to these data, turn the
            % option off so the Monte Carlo can run, if necessary. Also,
            % warn the user
            NullModelOpt = 'off';
            warning('The supplied null model does not apply to these data.')
        end
    end
    
    % If the user requested Monte Carlo, perform it
    if strcmp(MCOpt,'on') && strcmp(NullModelOpt,'off')
        
        % Figure out how many failures we're allowed
        nFailsThresh = ceil(MCpThresh * MCnSamples);
        iSample = 1;
        iFails = zeros([4,1]);
        MCInfoVals = NaN([4,MCnSamples]);
        ny = length(unique(y));
        nx1 = length(unique(x1));
        nx2 = length(unique(x2));
        ly = length(y);
        
        % Perform the Monte Carlo Trials
        while (nnz(iFails < nFailsThresh) > 0) && (iSample <= MCnSamples)
            MCInfoVals(:,iSample) = PID(accumarray({y,x1(randperm(ly)),x2(randperm(ly))},ones([ly,1]),[ny,nx1,nx2]),SourcesMat,SetsMat,TransMat);
            iFails(MCInfoVals(:,iSample) >= InfoVal) = iFails(MCInfoVals(:,iSample) >= InfoVal) + 1;
            iSample = iSample + 1;
        end
        
        % Calculate the p-value
        p = iFails / (iSample - 1);
        
        % Correct for the resolution of the Monte Carlo trials
        p(p == 0) = 1/(2*MCnSamples);
        
        % Record that the Monte Carlo method was used
        s = 0;
    end
    
    
elseif strcmp(Method,'3PID')
    
    %% Perform the 3-Variable PID Calculation
    
    % Parse the inputs into trial states
    if STFlag
        t1 = cell2mat(VariableIDs(:,3)) - min(cell2mat(VariableIDs(:,3)));
        t2 = cell2mat(VariableIDs(:,3)) - max(cell2mat(VariableIDs(:,3)));
        y = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},(t1(1) + 1):(end + t2(1))));
        x1 = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},(t1(2) + 1):(end + t2(2))));
        x2 = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},(t1(3) + 1):(end + t2(3))));
        x3 = squeeze(StatesRaster{VariableIDs{4,1}}(VariableIDs{4,2},(t1(4) + 1):(end + t2(4))));
    else
        y = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},VariableIDs{1,3},:));
        x1 = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},VariableIDs{2,3},:));
        x2 = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},VariableIDs{3,3},:));
        x3 = squeeze(StatesRaster{VariableIDs{4,1}}(VariableIDs{4,2},VariableIDs{4,3},:));
    end
    
    % Error check the integer states
    if ~isequal(unique(x1),(1:length(unique(x1))))
        [B,I,x1] = unique(x1);
    end
    if ~isequal(unique(x2),(1:length(unique(x2))))
        [B,I,x2] = unique(x2);
    end
    if ~isequal(unique(x3),(1:length(unique(x3))))
        [B,I,x3] = unique(x3);
    end
    if ~isequal(unique(y),(1:length(unique(y))))
        [B,I,y] = unique(y);
    end
    
    % Make a counts matrix
    Counts = accumarray({y,x1,x2,x3},ones(size(y)),[length(unique(y)),length(unique(x1)),length(unique(x2)),length(unique(x3))]);
    
    % Load the necessary matrices
    load('3PIDMats.mat')
    
    % Calculate the pairwise mutual information
    InfoVal = PID(Counts,SourcesMat,SetsMat,TransMat);
    
    % If the user requested Monte Carlo, perform it
    if strcmp(MCOpt,'on')
        
        % Figure out how many failures we're allowed
        nFailsThresh = ceil(MCpThresh * MCnSamples);
        iSample = 1;
        iFails = zeros([18,1]);
        MCInfoVals = NaN([18,MCnSamples]);
        ny = length(unique(y));
        nx1 = length(unique(x1));
        nx2 = length(unique(x2));
        nx3 = length(unique(x3));
        ly = length(y);
        
        % Perform the Monte Carlo Trials
        while (nnz(iFails < nFailsThresh) > 0) && (iSample <= MCnSamples)
            MCInfoVals(:,iSample) = PID(accumarray({y,x1(randperm(ly)),x2(randperm(ly)),x3(randperm(ly))},ones([ly,1]),[ny,nx1,nx2,nx3]),SourcesMat,SetsMat,TransMat);
            iFails(MCInfoVals(:,iSample) >= InfoVal) = iFails(MCInfoVals(:,iSample) >= InfoVal) + 1;
            iSample = iSample + 1;
        end
        
        % Calculate the p-value
        p = iFails / (iSample - 1);
        
        % Correct for the resolution of the Monte Carlo trials
        p(p == 0) = 1/(2*MCnSamples);
        
        % Record that the Monte Carlo method was used
        s = 0;
    end
    
    
elseif strcmp(Method,'3TE')
    
    %% Perform the 3-Variable Transfer Entropy Calculation
    
    % Error check the times
    if VariableIDs{2,3} >= VariableIDs{1,3}
        error('Past state of receiver after the future state of the receiver!')
    elseif VariableIDs{3,3} >= VariableIDs{1,3}
        error('Past state of transmitter 1 after the future state of the receiver!')
    elseif VariableIDs{4,3} >= VariableIDs{1,3}
        error('Past state of transmitter 2 after the future state of the receiver!')
	elseif ~isequal(cell2mat(VariableIDs(1,1:2)),cell2mat(VariableIDs(2,1:2)))
        error('Past and future states of the receiver come from different variables!')
    end
    
    % Parse the inputs into trial states
    if STFlag
        t1 = cell2mat(VariableIDs(:,3)) - min(cell2mat(VariableIDs(:,3)));
        t2 = cell2mat(VariableIDs(:,3)) - max(cell2mat(VariableIDs(:,3)));
        y = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},(t1(1) + 1):(end + t2(1))));
        x1 = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},(t1(2) + 1):(end + t2(2))));
        x2 = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},(t1(3) + 1):(end + t2(3))));
        x3 = squeeze(StatesRaster{VariableIDs{4,1}}(VariableIDs{4,2},(t1(4) + 1):(end + t2(4))));
    else
        y = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},VariableIDs{1,3},:));
        x1 = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},VariableIDs{2,3},:));
        x2 = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},VariableIDs{3,3},:));
        x3 = squeeze(StatesRaster{VariableIDs{4,1}}(VariableIDs{4,2},VariableIDs{4,3},:));
    end
    
    % Error check the integer states
    if ~isequal(unique(x1),(1:length(unique(x1))))
        [B,I,x1] = unique(x1);
    end
    if ~isequal(unique(x2),(1:length(unique(x2))))
        [B,I,x2] = unique(x2);
    end
    if ~isequal(unique(x3),(1:length(unique(x3))))
        [B,I,x3] = unique(x3);
    end
    if ~isequal(unique(y),(1:length(unique(y))))
        [B,I,y] = unique(y);
    end
    
    % Make a counts matrix
    Counts = accumarray({y,x1,x2,x3},ones(size(y)),[length(unique(y)),length(unique(x1)),length(unique(x2)),length(unique(x3))]);
    
    % Load the necessary matrices
    load('3PIDMats.mat')
    load('TE3Redux.mat')
    
    % Calculate the pairwise mutual information
    InfoVal = TE3Redux*PID(Counts,SourcesMat,SetsMat,TransMat);
    
    % If the user requested Monte Carlo, perform it
    if strcmp(MCOpt,'on')
        
        % Figure out how many failures we're allowed
        nFailsThresh = ceil(MCpThresh * MCnSamples);
        iSample = 1;
        iFails = zeros([4,1]);
        MCInfoVals = NaN([4,MCnSamples]);
        ny = length(unique(y));
        nx1 = length(unique(x1));
        nx2 = length(unique(x2));
        nx3 = length(unique(x3));
        ly = length(y);
        
        % Perform the Monte Carlo Trials
        while (nnz(iFails < nFailsThresh) > 0) && (iSample <= MCnSamples)
            MCInfoVals(:,iSample) = TE3Redux*...
                PID(accumarray({y,x1(randperm(ly)),x2(randperm(ly)),x3(randperm(ly))},ones([ly,1]),[ny,nx1,nx2,nx3]),SourcesMat,SetsMat,TransMat);
            iFails(MCInfoVals(:,iSample) >= InfoVal) = iFails(MCInfoVals(:,iSample) >= InfoVal) + 1;
            iSample = iSample + 1;
        end
        
        % Calculate the p-value
        p = iFails / (iSample - 1);
        
        % Correct for the resolution of the Monte Carlo trials
        p(p == 0) = 1/(2*MCnSamples);
        
        % Record that the Monte Carlo method was used
        s = 0;
    end
    
    
elseif strcmp(Method,'InfoGain')
    
    %% Perform the Information Gain Calculation
    
    % Error check the times
    if VariableIDs{2,3} >= VariableIDs{1,3}
        error('Past state of receiver after the future state of the receiver!')
	elseif ~isequal(cell2mat(VariableIDs(1,1:2)),cell2mat(VariableIDs(2,1:2)))
        error('Past and future states of the receiver come from different variables!')
    end
    
    % Parse the inputs into trial states
    if STFlag
        t1 = cell2mat(VariableIDs(:,3)) - min(cell2mat(VariableIDs(:,3)));
        t2 = cell2mat(VariableIDs(:,3)) - max(cell2mat(VariableIDs(:,3)));
        y = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},(t1(1) + 1):(end + t2(1))));
        x1 = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},(t1(2) + 1):(end + t2(2))));
        x2 = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},(t1(3) + 1):(end + t2(3))));
    else
        y = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},VariableIDs{1,3},:));
        x1 = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},VariableIDs{2,3},:));
        x2 = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},VariableIDs{3,3},:));
    end
    
    % Error check the integer states
    if ~isequal(unique(x1),(1:length(unique(x1))))
        [B,I,x1] = unique(x1);
    end
    if ~isequal(unique(x2),(1:length(unique(x2))))
        [B,I,x2] = unique(x2);
    end
    if ~isequal(unique(y),(1:length(unique(y))))
        [B,I,y] = unique(y);
    end
    
    % Make a counts matrix
    Counts = accumarray({x2,y,x1},ones(size(y)),[length(unique(x2)),length(unique(y)),length(unique(x1))]);
    
    % Load the necessary matrices
    load('2PIDMats.mat')
    
    % Calculate the pairwise mutual information
    PIDTerms = PID(Counts,SourcesMat,SetsMat,TransMat);
    
    InfoVal = MutualInfo(squeeze(sum(Counts,3))) - PIDTerms(1);
    
    % If the user requested Monte Carlo, perform it
    if strcmp(MCOpt,'on')
        
        % Figure out how many failures we're allowed
        nFailsThresh = ceil(MCpThresh * MCnSamples);
        iSample = 1;
        iFails = 0;
        MCInfoVals = NaN([1,MCnSamples]);
        ny = length(unique(y));
        nx1 = length(unique(x1));
        nx2 = length(unique(x2));
        ly = length(y);
        
        % Perform the Monte Carlo Trials
        while (nnz(iFails < nFailsThresh) > 0) && (iSample <= MCnSamples)
            PIDTerms = PID(accumarray({x2(randperm(ly)),y,x1},ones([ly,1]),[nx2,ny,nx1]),SourcesMat,SetsMat,TransMat);
            MCInfoVals(iSample) = MutualInfo(squeeze(sum(Counts,3))) - PIDTerms(1);
            iFails(MCInfoVals(iSample) >= InfoVal) = iFails(MCInfoVals(iSample) >= InfoVal) + 1;
            iSample = iSample + 1;
        end
        
        % Calculate the p-value
        p = iFails / (iSample - 1);
        
        % Correct for the resolution of the Monte Carlo trials
        p(p == 0) = 1/(2*MCnSamples);
        
        % Record that the Monte Carlo method was used
        s = 0;
    end
    
    
elseif strcmp(Method,'InfoTrans')
    
    %% Perform the Information Transfer
    
    % Error check the times
    if VariableIDs{2,3} >= VariableIDs{1,3}
        error('Past state of receiver after the future state of the receiver!')
    elseif VariableIDs{4,3} >= VariableIDs{1,3}
        error('Past state of transmitter after the future state of the receiver!')
	elseif ~isequal(cell2mat(VariableIDs(1,1:2)),cell2mat(VariableIDs(2,1:2)))
        error('Past and future states of the receiver come from different variables!')
    end
    
    % Parse the inputs into trial states
    if STFlag
        t1 = cell2mat(VariableIDs(:,3)) - min(cell2mat(VariableIDs(:,3)));
        t2 = cell2mat(VariableIDs(:,3)) - max(cell2mat(VariableIDs(:,3)));
        y = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},(t1(1) + 1):(end + t2(1))));
        x1 = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},(t1(2) + 1):(end + t2(2))));
        x2 = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},(t1(3) + 1):(end + t2(3))));
        x3 = squeeze(StatesRaster{VariableIDs{4,1}}(VariableIDs{4,2},(t1(4) + 1):(end + t2(4))));
    else
        y = squeeze(StatesRaster{VariableIDs{1,1}}(VariableIDs{1,2},VariableIDs{1,3},:));
        x1 = squeeze(StatesRaster{VariableIDs{2,1}}(VariableIDs{2,2},VariableIDs{2,3},:));
        x2 = squeeze(StatesRaster{VariableIDs{3,1}}(VariableIDs{3,2},VariableIDs{3,3},:));
        x3 = squeeze(StatesRaster{VariableIDs{4,1}}(VariableIDs{4,2},VariableIDs{4,3},:));
    end
    [waste1,waste2,x4] = unique(sub2ind([length(unique(x1)),length(unique(x3))],x1,x3));
    
    % Error check the integer states
    if ~isequal(unique(x1),(1:length(unique(x1))))
        [B,I,x1] = unique(x1);
    end
    if ~isequal(unique(x2),(1:length(unique(x2))))
        [B,I,x2] = unique(x2);
    end
    if ~isequal(unique(x4),(1:length(unique(x4))))
        [B,I,x4] = unique(x4);
    end
    if ~isequal(unique(y),(1:length(unique(y))))
        [B,I,y] = unique(y);
    end
    
    % Make a counts matrix
    Counts1 = accumarray({x2,y,x4},ones(size(y)),[length(unique(x2)),length(unique(y)),length(unique(x4))]);
    Counts2 = accumarray({x2,y,x1},ones(size(y)),[length(unique(x2)),length(unique(y)),length(unique(x1))]);
    
    % Load the necessary matrices
    load('2PIDMats.mat')
    
    % Calculate the pairwise mutual information
    PIDTerms1 = PID(Counts1,SourcesMat,SetsMat,TransMat);
    PIDTerms2 = PID(Counts2,SourcesMat,SetsMat,TransMat);
    
    InfoVal = PIDTerms1(1) - PIDTerms2(1);
    
    % If the user requested Monte Carlo, perform it
    if strcmp(MCOpt,'on')
        
        % Figure out how many failures we're allowed
        nFailsThresh = ceil(MCpThresh * MCnSamples);
        iSample = 1;
        iFails = 0;
        MCInfoVals = NaN([1,MCnSamples]);
        ny = length(unique(y));
        nx1 = length(unique(x1));
        nx2 = length(unique(x2));
        nx4 = length(unique(x4));
        ly = length(y);
        
        % Perform the Monte Carlo Trials
        while (nnz(iFails < nFailsThresh) > 0) && (iSample <= MCnSamples)
            RandX2 = x2(randperm(ly));
            PIDTerms1 = PID(accumarray({RandX2,y,x4},ones(size(y)),[nx2,ny,nx4]),SourcesMat,SetsMat,TransMat);
            PIDTerms2 = PID(accumarray({RandX2,y,x1},ones(size(y)),[nx2,ny,nx1]),SourcesMat,SetsMat,TransMat);
            MCInfoVals(iSample) = PIDTerms1(1) - PIDTerms2(1);
            iFails(MCInfoVals(iSample) >= InfoVal) = iFails(MCInfoVals(iSample) >= InfoVal) + 1;
            iSample = iSample + 1;
        end
        
        % Calculate the p-value
        p = iFails / (iSample - 1);
        
        % Correct for the resolution of the Monte Carlo trials
        p(p == 0) = 1/(2*MCnSamples);
        
        % Record that the Monte Carlo method was used
        s = 0;
    end
    
else
    error('An invalid information measurement method input was supplied.')
end
end


% Null Model Tester and Application
function p = nullpmf(InfoVal,Counts,nullModel,Method)
    
    % Check that the modes match
    if ~strcmp(nullModel.Method,Method)
        error('The method for the null model does not match the current method.')
    end
    
    % Check the number of variables
    if ndims(Counts) ~= length(nullModel.nCounts)
        error('The null model and these data do not have the same number of variables.')
    end
    
    % Check that the info values and pmf are cells
    if ~iscell(nullModel.info)
        nullModel.info = {nullModel.info};
    end
    if ~iscell(nullModel.pmf)
        nullModel.pmf = {nullModel.pmf};
    end
    
    % Check each variable to make sure that the number of states match
    for iVar = 1:ndims(Counts)
        tempCounts = Counts;
        for jVar = setdiff(1:ndims(Counts),iVar)
            tempCounts = sum(tempCounts,jVar);
        end
        if ~isequal(reshape(sort(squeeze(tempCounts),'ascend'),[length(tempCounts),1]),nullModel.nCounts{iVar})
            error(['The data and the null model do not have the same number of counts in variable ',num2str(iVar)])
        end
    end
    
    % Round the information values to correct for small rounding errors
    % associated with the system resolution (eps).
    InfoVal = round(InfoVal./nullModel.rFact)*nullModel.rFact;
    
    % Go through each information value and calculate the p-value
    p = NaN(size(InfoVal));
    for iMeas = 1:length(InfoVal)
        p(iMeas) = sum(nullModel.pmf{iMeas}(nullModel.info{iMeas} >= InfoVal(iMeas))) * nullModel.res;
    end
    
    % Correct for the resolution of the null model
    p(p == 0) = nullModel.res/2;
    
end
