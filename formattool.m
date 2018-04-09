%% FORMATTOOL - converts simple data format to dataraster format
% Organizes data from simple formats (e.g. voltage value vectors) to the
% dataraster format for processing in the information analysis. For each
% point in timelocks, data will be binned around the time lock point. One
% bin will start for time values equal to or greater than the timelock
% point and extending binsize time out in the future. Other bins will be
% added before and after the timelock point using identical boundary
% patterns. 
%
% Syntax: [formdata,timeboundaries] = formattool(data, time, timelocks, binsize, maxLead, maxLag, method)
%
% Input:
%   data (1 by nObs double array) - basic data to be converted. The data
%       can be continuous data (e.g. voltage values), or it can be discrete
%       (e.g. integer stimulus or behavior states, integer spike states (1
%       = spike). nObs is the total number of observations.
%   time (1 by nObs double array) - time values corresponding to the basic
%       data values in the data array. nObs is the total number of
%       observations. All variables using time must have the same units.
%   timelocks (1 by nTrials double array) - time lock points. These times will be used
%       allign sections of data in the dataraster. nTrials is the total number of
%       trials. All variables using time must have the same units.
%   binsize (double) - the size of the time bins. All variables using time
%       must have the same units.
%   maxLead (double) - the maximum amount of time analyzed prior to the
%       time lock points. All variables using time must have the same
%       units.
%   maxLag (double) - the maximum amount of time analyzed after the time
%       lock points. All variables using time must have the same units.
%   method (string) - the method to be used to bin the data. The method is 
%       only relevant if the bin size is changed. Possible choices are:
%           'count': total count of the data values in a bin. This is most
%               useful for spike counts where spikes are noted as 1 in the
%               data array.
%           'ave': average data values. This can be useful for continuous
%               valued data (e.g. voltage data).
%
% Outputs:
%   formdata (1 by nbins by nTrials double array) - the original data in
%       the dataraster format. The data can be easily added to a larger 
%       array with other data to form a large dataraster, which can then be 
%       processed by data2states.
%   timeboundaries (nbins by 2 double array) - the boundaries for the
%       corresponding time bin in formdata. The first column represents the
%       earliest time in the bin and the second column represents the
%       latest time in the bin. For instance, if timeboundaries(5,:) =
%       [1,2], formdata(1,5,8) includes data at time points greater than or 
%       equal to 1 time units, but less than 2 time units after 
%       timelocks(8).
%
% Examples:
%   Noisy Sin-Wave Data
%       time = 1:20000;
%       data = sin(2*pi*time/500) + 0.2*randn(size(time));
%       timelocks = 125:500:20000;
%       [formdata] = formattool(data, time, timelocks, 25, 100, 400, 'ave');
%   Sin-Wave Driven Spike Data
%       time = 1:20000;
%       x = 0.1*(sin(2*pi*time/500) + 1);
%       spikes = (x > rand(size(x)));
%       time(~spikes) = [];
%       data = ones(size(time));
%       timelocks = 125:500:20000;
%       [formdata] = formattool(data, time, timelocks, 25, 100, 400, 'count');
%   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Nicholas Timme
% Email: nicholas.m.timme@gmail.com
% February 2016; Last revision: 12-May-2016


function [formdata,timeboundaries] = formattool(data, time, timelocks, binsize, maxLead, maxLag, method)
%% Organize the Data

% Error check the time and data vectors
if length(data) ~= length(time)
    error('data and time must be equal length')
end

% Figure out the bin boundaries
TimeBounds = (-floor(maxLead/binsize):floor(maxLag/binsize))*binsize;

% Preallocate the formatted data
nBins = length(TimeBounds) - 1;
nTrials = length(timelocks);
formdata = NaN([1,nBins,nTrials]);
timeboundaries = NaN([nBins,2]);

% Make the timeboundaries list
for iBin = 1:nBins
    timeboundaries(iBin,:) = [TimeBounds(iBin),TimeBounds(iBin + 1)];
end

% Process the data
if strcmp(method,'count')
    
    for iTrial = 1:nTrials
        for iBin = 1:nBins
            formdata(1,iBin,iTrial) = sum(data((time >= (timelocks(iTrial) + TimeBounds(iBin))) & (time < (timelocks(iTrial) + TimeBounds(iBin + 1)))));
        end
    end
    
elseif strcmp(method,'ave')
    
    for iTrial = 1:nTrials
        for iBin = 1:nBins
            formdata(1,iBin,iTrial) = mean(data((time >= (timelocks(iTrial) + TimeBounds(iBin))) & (time < (timelocks(iTrial) + TimeBounds(iBin + 1)))));
        end
    end
    
end





end