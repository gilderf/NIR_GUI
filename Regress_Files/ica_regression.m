function [ r2train,r2test,rmsecvC,rmsepC  ] = ica_regression( x_train, y_train, x_test , y_test , range)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
addpath FastICA_25


% %% This calculates the ICA. We could tell the function how many ICA's we
% % want to calculate. If we don't, it will estimate how many real components
% % are found in the data set.
% [ica, A, W] = fastica(x_train);
% 
% % ica = the component spectra
% % W = weights
% 
% %% Let's plot the first ICA "spectrum".
% % figure, plot(range,ica(1,:))
% % xlabel('Wavelength')
% % title('ICA1 COEFF')
% 
% %% Let's plot the second ICA "spectrum".
% % figure, plot(range,ica(2,:))
% % xlabel('Wavelength')
% % title('ICA2 COEFF')
% 
% %% Correlations of individual ICA PCs and instant modulus
% 
% corr(A(:,1),y_train)
% corr(A(:,2),y_train)
% corr(A(:,3),y_train)
% corr(A(:,4),y_train)
% corr(A(:,5),y_train)
% 
% % not very good correlations individually here, either.

%% my try at ICA regression (based on PCR example). 
% seems to give worse results than PCR, I don't know if this is correct. In
% any case, this is just a first trial.
X = x_train;
y = y_train;
[n,p] = size(X);

[ica, A, W] = fastica(x_train);
% in ICA regression, we combine the score values of multiple components.
% for example, the components 1-6
% betaICA = regress(y-mean(y), W(1:10,:)');

betaICA = regress(y-mean(y), A(:,1:10));
betaICA = ica(1:10,:)'*betaICA;
betaICA = [mean(y) - mean(X)*betaICA; betaICA];

% betaICA = regress(y-mean(y), A);
% betaICA = ica*betaICA;
% betaICA = [mean(y) - mean(X)*betaICA; betaICA];


yfitICA_train = [ones(n,1) X]*betaICA;
yfitICA_test = [ones(n,1) x_test]*betaICA;

corr(yfitICA,y)



r2train=0;
r2test=0;
rmsecvC=0;
rmsepC=0

ncompmax= 10

for ncomp=1:ncompmax

betaICA = regress(y-mean(y), A(:,1:ncomp));
betaICA = ica(1:10,:)'*betaICA;
betaICA = [mean(y) - mean(X)*betaICA; betaICA];

yfitICA_train = [ones(n,1) X]*betaICA;
yfitICA_test = [ones(n,1) x_test]*betaICA;


s1 = sum((y_train - yfitICA_train).^2);
s2 = sum((y_train - mean(y_train)).^2);
r2train(ncomp) = (1 - (s1/s2))*100;
rmsepC(ncomp) = sqrt(mean((stats.Yresiduals).^2));
rmsecvC(ncomp) = sqrt(mean((y_train-yfitICA_train).^2));

s11 = sum((y_test - yfitICA_test).^2);
s22 = sum((y_test - mean(y_test)).^2);
r2test(ncomp) = (1 - (s11/s22))*100;
end

% Print the results Grid

figure;
title('Baseline per manuscript');
subplot(2,2,1);
plot(1:ncompmax,r2train)
xlabel('Lambda Value');
ylabel('R^2 Train');
hold on
plot(find(r2train == min(r2train)),min(r2train),'o')
plot(find(r2test == max(r2test)),r2train,'X')
hold off

subplot(2,2,3);
plot(1:ncompmax,r2test)
xlabel('Lambda Value');
ylabel('R^2 Test');
hold on
plot(find(r2test == max(r2test)),max(r2test),'o')
hold off

subplot(2,2,4);
plot(1:ncompmax,rmsepC)
xlabel('Lambda Value');
ylabel('rmsepC');
hold on
plot(find(rmsepC == min(rmsepC)),min(rmsepC),'o')
plot(find(r2test == max(r2test)),rmsepC,'X')
hold off

subplot(2,2,2);
plot(1:ncompmax,rmsecvC)
xlabel('Lambda Value');
ylabel('RMSECV');
hold on
plot(find(rmsecvC == min(rmsecvC)),min(rmsecvC),'o')
plot(find(r2test == max(r2test)),rmsecvC,'X')
hold off
title_to_display = strcat('Ridge Regression : ',func_param);
suptitle(title_to_display)

end

