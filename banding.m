function [ b ] = banding( p,k )
%UNTITLED3 Summary of this function goes here
% Generate a banding operator with given dimention and tuning parameter. Multiplying 
% it on a covariance matrix by componentwise product can provide a regularized estimator
% with the banding method.
% xamples{
% p <- 5;
% W <- banding(p,k=2) ;
% W;
% Bickel, P and Levina, E, Regularized estimation of large covariance matrices,
% Annals of Statistics, 36, 199-227 (2008).

k = min(k,2*(p-1));
%w = matrix(0,p,p);
for i = 1:p
    
    for j = 1:p
        
            b(i,j) = 1 * (abs(i-j) <= k);
            
    end
    
end

end

