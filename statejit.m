%% statejit - jitter states
% Jitters states with adjacent states in time order, trial order, or
% whatever order they are presented. The jittering is performed using a
% Gaussian distribution with a supplied standard deviation (in units of
% bins). 
%
% Syntax: [randstates] = statejit(states,std)
%
% Input:
%   states (vector double) - 
%   std (scalar) - 
%
%
% Outputs:
%   randstates (vector double) - randomized version of the input states
%
%
% Examples:
% Jitter with 2 bin standard deviation
% states = rand([100,1]);
% randstates = statejit(states,2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% August 2016; Last revision: 01-Aug-2016


function [randstates] = statejit(states,stdev)

[waste,I] = sort((1:length(states)) + stdev*randn([1,length(states)]));
randstates = states(I);
