function [TermLabels] = PIDTermLabels(N)
%PIDTermLabels creates a list of the term labels to ease interpretation of
%the PID output.
%   [TermLabels] = PIDTermLabels(N) is a cellular array that contains the
%   term labels for the output of the PID function. It follows the notation
%   of Williams and Beer.
%
%   P. L. Williams and R. D. Beer, arXiv:1004.2515v1 (2010).
%
%   Inputs
%
%   N: The number of X variables in the PID calculation.
%
%   Outputs
%
%   TermLabels: A cellular array that contains the labels for the PID
%   terms.
%
%
%       Version 2.0

% Version Information
%
%   1.0: 10/6/11 - Original program created before and modified up to this
%   date. (Nick Timme)
%
%   2.0: 3/20/13 - Program formatting modified for inclusion in the
%   toolbox. Also, printing term labels was replaced by a cellular array.
%   (Nick Timme)
%

%==============================================================================
% Copyright (c) 2013, The Trustees of Indiana University
% All rights reserved.
% 
% Authors: Nick Timme (nmtimme@umail.iu.edu)
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
%   1. Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
% 
%   2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
% 
%   3. Neither the name of Indiana University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==========================================================================


% Obtain necessary lattice arrays from Lattice.
[SourcesMat,SetsMat,TransMat] = PIDLattice(N);

% First, label all of the sources. 

[NumSets,NumSources]=size(SetsMat);

for i=1:NumSources
    Var=find(SourcesMat(i,:)==1);
    templabel=['{'];
    for j=1:length(Var)
        templabel=[templabel,num2str(Var(j))];
    end
    templabel=[templabel,'}'];
    eval(['Source',num2str(i),'Label=templabel;']);
end

% Now label each set
TermLabels = cell([NumSets,1]);
for i=1:NumSets
    Sources=find(SetsMat(i,:)==1);
    templabel=[''];
    for j=1:length(Sources)
        eval(['templabel=[templabel,Source',num2str(Sources(j)),'Label];']);
    end
    TermLabels{i} = templabel;
end




end

