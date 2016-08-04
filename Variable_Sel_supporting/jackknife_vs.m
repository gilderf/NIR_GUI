function [jackknife_variables, p] = jackknife_vs(X,Y,ncomp,threshold,perturbation_parameter)

if nargin < 5
    perturbation_parameter = 0;
end

if nargin < 4
    threshold = 0.05;
end

KFold = numel(Y);
Indices = crossvalind('Kfold', numel(Y), KFold);

yfit_PLS = zeros(size(Y));
beta_cross_val = zeros(numel(Y),size(X,2)+1);

for ii=1:KFold
    ind_val = find(Indices == ii);
    ind_pred = find(Indices ~= ii);

    [XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(X(ind_pred,:),Y(ind_pred),ncomp);
    beta_cross_val(ii,:) = betaPLS;
    yfit_PLS(ind_val) = [ones(numel(ind_val),1) X(ind_val,:)]*betaPLS;
end

n = numel(Y);
betaHat = (1/n)*sum(beta_cross_val,1);
varBeta = ((n - 1)/ n).*sum((beta_cross_val - repmat(betaHat,[n 1])).^2);
% T = betaHat./sqrt(varBeta);
T = betaHat./(sqrt(varBeta) + perturbation_parameter);

p = tcdf(T,n);

jackknife_variables = find(p < threshold);