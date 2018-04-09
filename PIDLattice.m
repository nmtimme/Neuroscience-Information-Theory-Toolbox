function [SourcesMat,SetsMat,TransMat] = PIDLattice(N)
%PIDLattice creates necessary lattice variables for PID calculations.
%   [SourcesMat,SetsMat,TransMat] = PIDLattice(N) is a collection of
%   variables that describe the partial information lattice. These
%   variables depend on the number of X variables that are to be used in
%   the PID calculation. The purpose of this program is to speed up the
%   process of performing numerous PID calculations with the same number of
%   X variables. Lattice must be run prior to using PID. It is based on 
%   work by Williams and Beer.
%
%   P. L. Williams and R. D. Beer, arXiv:1004.2515v1 (2010).
%
%   N is the number of X variables that will be considered.
%
%   Inputs
%
%   N: The number of X variables that will be considered in the PID
%   calculation. 
%
%   Outputs
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
%
%       Version 2.0

% Version Information
%
%   1.0: 12/6/11 - Original program created and modified up to this date.
%   (Nick Timme)
%
%   2.0: 3/20/11 - Formatting of program modified for inclusion in the new
%   toolbox. (Nick Timme)

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




% The number of possible sources is simply 2^N-1 (we'll cut the empty set
% out in a minute)


SourcesMat=zeros([(2^N),N]);
Source=zeros([1,N]);
for i=1:(2^N)
    for j=1:N
        Source(j)=rem(ceil(i/(2^(j-1))),2);
    end
    SourcesMat(i,:)=Source;
end


% Now cut out the empty set.

SourcesMat(2^N,:)=[];


% Now we go through the possible sets of sources to find the correct ones

% Set the maximum number of sets of sources. The program will go past this
% if necessary, but this pre-allocates some space. For N=4, there are 
% around 150 sets.
MaxSize=300; 

% pre-allocate space for the largest possible number of sets
SetsMat=zeros([MaxSize,(2^N-1)]); 

% pre-allocate space for a possible set
PossibleSet=zeros([1,(2^N-1)]); % pre-allocate space for a possible set

% create a running count of the number of sets of sources we've found
index=1; 

for i=1:2^(2^N-1)
    for j=1:(2^N-1) % create the possible set of sources
        PossibleSet(j)=rem(ceil(i/(2^(j-1))),2);
    end
    % create an indicator that will be tripped if the set of sources fails
    % to meet some requirement
    Indicator=0; 
    if sum(PossibleSet)>1 % sets of sources with only one source are all legit
        temp=find(PossibleSet==1);
        for j=1:(length(temp)-1)
            for k=(j+1):length(temp)
                if (dot(SourcesMat(temp(j),:),SourcesMat(temp(k),:))...
                        ==sum(SourcesMat(temp(j),:)))||...
                        (dot(SourcesMat(temp(j),:),...
                        SourcesMat(temp(k),:))==sum(SourcesMat(temp(k),:)))
                    Indicator=1;
                end
            end
        end
    elseif sum(PossibleSet)==0 % The empty set is not legit
        Indicator=1;
    end
    if Indicator==0 % If the indicator hasn't been tripped, the set is good
        SetsMat(index,:)=PossibleSet;
        index=index+1;
    end
end

% Eliminate the extra space in the matrix of sets
SetsMat(index:MaxSize,:)=[];


[r,waste]=size(SetsMat);

% Now we must make a matrix that contains information about the connections
% between the sets on the lattice

ConMat=zeros([r,r]);


for i=1:r
    for j=1:r
        if ~isequal(i,j)
            Indicator=0;
            temp1=find(SetsMat(i,:)==1); %list of sources in first set
            temp2=find(SetsMat(j,:)==1); %list of sources in second set
            for k=1:length(temp1)
                Indicator1=0;
                for l=1:length(temp2)
                    temp3=dot(SourcesMat(temp1(k),:),...
                        SourcesMat(temp2(l),:));
                    if temp3==sum(SourcesMat(temp2(l),:))
                        Indicator1=1;
                    end
                end
                if Indicator1==1
                    Indicator=Indicator+1;
                end
            end
            if Indicator==length(temp1)
                ConMat(i,j)=1;
            end
        end
    end
end



% Now construct the transform matrix

TransMat=zeros([r,r]);

[NumSets,waste]=size(SetsMat);

% Create a list of found states. We will need those so we can proceed "up"
% the lattice.
Found=zeros([NumSets,1]); 
FinishedFound=ones([NumSets,1]);

% First, find the bottom node. This node will not have any connections
% below it on the lattice.

temp=ConMat*ones([NumSets,1]);

Found(temp==0)=1;

% The partial information term for this node will be equal to its minimum
% information.
TransMat(Found==1,Found==1)=1;

% Now go through the other nodes on the lattice
NodesNeeded=sum(ConMat,2); 

while ~isequal(Found,FinishedFound)
    temp1=zeros([NumSets,1]);
    temp1(NodesNeeded<=sum(Found))=1;
    temp1=temp1-Found;
    NewFinds=find(temp1==1);
    for i=1:length(NewFinds)
        % The found list might have more things than we need
        if dot(ConMat(NewFinds(i),:)',Found)==sum(ConMat(NewFinds(i),:)) 
            
            % The output will always need the minimum information value for
            % that given set of sources
            TransMat(NewFinds(i),NewFinds(i))=1; 
            
            % Find all the partial information terms that are lower on the
            % matrix that must be subtracted
            ToSubtract=find(ConMat(NewFinds(i),:)==1); 
            for j=1:length(ToSubtract)
                TransMat(NewFinds(i),:)=TransMat(NewFinds(i),:)-TransMat(ToSubtract(j),:);
            end
            Found(NewFinds(i))=1;
        end
    end
end



end

