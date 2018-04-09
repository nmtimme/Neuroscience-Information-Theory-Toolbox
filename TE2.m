function [TE] = TE2(CountsMat)
% [TE] = TE2(CountsMat) TE calculates the transfer entropy from the X2
% variable to the Y variable assuming the X1 variable is the history of the
% signal or process that preceeded Y. It is based on the transfer entropy
% paper by Schrieber:
%
% T. Schreiber, Measuring Information Transfer, PRL 85 (2000).
%
% CountsMat is a tensor that contains the counts of the various states of the
% variables. The first index corresponds to the state of the Y variable.
% The second through N+1 indexes correspond to the states of the X1 to XN
% variables. New to version 2 is the ability to use counts matrices with an
% N+2 index that corresponds to a different data set.
%
% This program was written by Nicholas Timme. 
% 
% Version 1 was last updated on 5/15/12. 
% 
% Version 2 was written on 8/22/12. It incorporates the ability to
% process many identically sized counts matrices simultaneously via an
% extra dimension on the counts matrix.

[ny,nx1,nx2,nd]=size(CountsMat);

TotalCounts=sum(sum(sum(CountsMat,1),2),3);

P=CountsMat./TotalCounts(ones([1,ny]),ones([1,nx1]),ones([1,nx2]),:);

Px1x2=sum(P,1);

PyGivx1x2=P./Px1x2(ones([1,ny]),:,:,:);

Px1=sum(sum(P,1),3);
Pyx1=sum(P,3);

PyGivx1=Pyx1(:,:,ones([1,nx2]),:)./Px1(ones([1,ny]),:,ones([1,nx2]),:);


TE=P.*log2(PyGivx1x2./PyGivx1);
TE(~isfinite(TE))=0;
TE=squeeze(sum(sum(sum(TE,1),2),3));





end

