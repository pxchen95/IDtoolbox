function [mledim] = mledim(X,k1,k2)
% MLE estimate of intrinsic dimension
% The inputs are:
% X = d x N data matrix, (each column = 1 data pt)
% where d is the dimension and N is the number of observations
% k1 = first k (# nearest neighbors) for which to compute estimator
% k2 = last k for which to compute estimator
% The final estimate is averaged over k = k1...k2
%
% Default values:  
% mledim(X,k)=mledim(X,k,k)
% mledim(X) = mledim(X,10)
%
% Reference:  E. Levina and P.J. Bickel (2005).  
% "Maximum Likelihood Estimation  of Intrinsic Dimension."  
% In Advances in NIPS 17, Eds. L. K. Saul, Y. Weiss, L. Bottou. 
%
% Note: There is a small error in the paper which has been corrected 
% in this code: the normalizer should be k-2, not k-1. 

% check values of k1, k2
if(nargin==1) %k1, k2 not provided. Setting k1=k2=10.
    k1=10;           
    k2=10;           
end
if(nargin==2) %k2 not provided.  Setting k2=k1.
    k2=k1;           
end
if(k1<3) %k1 < 3.  Setting k1=10. 
     k1=10;          
end
if(k2<k1)  %k2 < k1.  Setting k2 = k1.
    k2=k1;           
end

% Compute matrix of log nearest neighbor distances
[~, N] = size(X);
X2 = sum(X.^2,1); 
knnmatrix= zeros(k2,N);

% To prevent memory problems on smaller machines, for large sample sizes 
% neighbors are found one point at a time

if (N < 4000)
    distance = sort(repmat(X2,N,1)+repmat(X2',1,N)-2*(X')*X);
    knnmatrix= .5*log(distance(2:k2+1,:));
else
    for i = 1:N
	      distance=sort(repmat(X2(i),1,N)+X2-2*X(:,i)'*X);
              knnmatrix(:,i)= .5*log(distance(2:k2+1))'; 
    end
end  

% Compute the estimate
S = cumsum(knnmatrix,1);
indexk = repmat((k1:k2)', 1, N);
dhat = -(indexk-2)./(S(k1:k2, :) - knnmatrix(k1:k2,:).*indexk);

% Take average over observations and over k1...k2                           
mledim = mean(mean(dhat));  