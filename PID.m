function [PITerms] = PID(CountsMat,SourcesMat,SetsMat,TransMat)
%PID calculates the partial information decomposition.
%   [PITerms] = PID(CountsMat,SourcesMat,SetsMat,TransMat) is an array that
%   contains the partial information terms. It requires several arrays that
%   are calculated using PIDLattice. This program is based on the work of 
%   Williams and Beer.
%
%   P. L. Williams and R. D. Beer, arXiv:1004.2515v1 (2010).
%
%   Inputs
%
%   CountsMat: An array that contains the counts (or joint probability 
%   values) of the various states of the variables. The first index 
%   corresponds to the state of the Y variable. The second through N+1 
%   indexes correspond to the states of the X1 to XN variables. 
%
%   SourcesMat: An array that contains all possible sources for a given set
%   of X variables. Each row is a source and each column is a variable. 
%
%   SetsMat: An array that contains the sets of sources in script A. Each 
%   row is a set and each column is a source. Sources are numbered 
%   according to SourcesMat. i.e. column 1 is source 1 (row 1 on 
%   SourcesMat), column 2 is source 2 (row 2 on SourcesMat), etc.
%
%   TransMat: An array that transforms the minimum informations directly 
%   to the minimum information terms. The dimensions are r by r (where r 
%   is the number of sets of sources). The partial information terms will 
%   be found via simple matrix multiplication TransMat*MinimumInfos.
%
%   Outputs
%
%   PITerms: An array that contains the partial information decomposition
%   term values in bits.
%
%
%       Version 2.0

% Version Information
%
%   1.0: 4/26/12 - The original version of the program was created before
%   and modified up to this data. (Nick Timme)
%
%   2.0: 3/20/13 - The program's formatting was modified for inclusion in
%   the toolbox. (Nick Timme)
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





% Find the number of states of Y and the X's
SizeCountsMat=size(CountsMat);

% Find the number of sources and sets of sources
[NumSets,NumSources]=size(SetsMat);

% Create a specific information matrix
SpecificInfos=zeros([SizeCountsMat(1),NumSources]);

% Reshape the counts matrix
    
NewCountsMat=reshape(CountsMat,SizeCountsMat(1),prod(SizeCountsMat)/...
    SizeCountsMat(1));
    
% Find the probabilities for each state of Y

TotalCounts=sum(sum(NewCountsMat));
YCounts=sum(NewCountsMat,2);
Py=YCounts/TotalCounts;

% Now we must go through each source and calculate the appropriate
% probabilities and the specific information

for i=1:NumSources
    TempCountsMat=CountsMat;
    
    % These are the R variables that aren't in the source
    ToElim=find(SourcesMat(i,:)==0); 
    
    % Sum over the variables that aren't in the source
    for j=1:length(ToElim)
        TempCountsMat=sum(TempCountsMat,ToElim(j)+1);
    end
    if size(TempCountsMat,1) ~= 1
        TempCountsMat=squeeze(TempCountsMat);
    else
        TempCountsMat=squeeze(TempCountsMat);
        TempCountsMat=reshape(TempCountsMat,[1,size(TempCountsMat)]);
    end
    % Now reshape the counts matrix to be Y states by A (source) states (A
    % states are combinations of X variables)
    TempSizeCountsMat=size(TempCountsMat);
    TempCountsMat=reshape(TempCountsMat,TempSizeCountsMat(1),...
        prod(TempSizeCountsMat)/TempSizeCountsMat(1));
    
    % Now we must calculate the conditional probabilities
    TempSizeCountsMat=size(TempCountsMat);
    
    % First, a conditional on y
    Pay=zeros(TempSizeCountsMat);
    for j=1:TempSizeCountsMat(1)
        Pay(j,:)=TempCountsMat(j,:)/YCounts(j);
    end
    
    % Second, y conditional on a
    Pya=zeros(TempSizeCountsMat);
    aCounts=sum(TempCountsMat,1);
    for j=1:TempSizeCountsMat(2)
        Pya(:,j)=TempCountsMat(:,j)/aCounts(j);
    end
    
    % Reset infinite values to zero
    Pay(~isfinite(Pay))=0;
    Pya(~isfinite(Pya))=0;
    
    % Calculate the specific informations
    
    for j=1:TempSizeCountsMat(1)
        temp1=log2(Pya(j,:)/Py(j));
        
        % If Pya/Py=0, the log will blow up, but in that case Pay=0, so the
        % limit will be 0. If Pya=0 and Py=0, we want the log to go to 0. 
        % Note, Pya>0 and Py=0 is not possible (or at least shouldn't be).
        temp1(~isfinite(temp1))=0; 
        SpecificInfos(j,i)=Pay(j,:)*temp1';
    end
end



% Create a matrix of minimum information values

MinInfos=zeros([NumSets,1]);

% Iterate through each set of sources.

for i=1:NumSets
    
    % select the relevant specific informations for each set of sources
    RelSpecInfos=SpecificInfos(:,SetsMat(i,:)==1); 
    
    % dot the minimum infos with the matching probability
    MinInfos(i)=dot(min(RelSpecInfos,[],2),Py); 
end


% Now calculate the partial information terms

PITerms=TransMat*MinInfos;


end

