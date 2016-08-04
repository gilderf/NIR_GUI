function [ w ] = tapering( p,k )
%UNTITLED2 Summary of this function goes here
% Generate a tapering operator with given dimention and tuning parameter. Multiplying 
% it on a covariance matrix by componentwise product can provide a regularized estimator
% with the tapering method.
% references{
% Cai, T, Zhang, CH and Zhou, H, Optimal rates of convergence for covariance
% matrix estimation, Annals of Statistics, 38, 2118-2144 (2010).

k = min(k,2*(p-1));
%w = matrix(0,p,p);
for i = 1:p
    
    for j = 1:p
        
            if (abs(i-j)<= k/2)
                
            w(i,j) = 1;
           
            end
            
            if ((k/2 < abs(i-j)) && (abs(i-j) < k))
                
                w(i,j) = 2 - 2* abs(i-j)/k;
                
            end
            
    end
    
end
%return(w);
end

