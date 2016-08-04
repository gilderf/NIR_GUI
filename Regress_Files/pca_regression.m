function [ r2train,r2test,rmsecvC,rmsepC  ] = pca_regression( x_train, y_train, x_test , y_test ,range, func_param);
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

[COEFF,SCORE,latent,tsquare] = princomp(x_train);
% COEFF contains the new variables calculated by PCA. We can consider them
% as PCA spectra and plot them to see if we can link them to certain
% constituents. PCA calculates as many new variables as we have in the
% original data set. However, in reality, only a few first PC:s are
% investigated, as the first PC explains the most of the variance in the data
% set, the second PC the second most on, and so on. So we need usually 1-10
% PC:s, as the rest contain only noise.

% SCORE contains the weights of each new variable for each original
% spectrum. So, if a PC would be directly linked to i.e. collagen, the
% score of this component would tell about the collagen content in the
% sample. In reality, we rarely find a single PC that could be used.

% latent tells us the variance accounted for by each component, so if we
% display the cumulative sum for the first 10 components, it reveals that
% including 10 component explains 73%
100*cumsum(latent(1:10))./sum(latent);

%% Let's plot the PC1 "spectrum".
% Savitzky-Golay smoothing shortens the spectra by two variables,
% therefore, the wavelength axis is modified accordingly.
% figure, plot(range,COEFF(1,:))
% xlabel('Wavelength')
% title('PC1 COEFF')
%
% %% Let's plot the PC2 "spectrum".
% figure, plot(range,COEFF(:,1))
% xlabel('Wavelength')
% title('PC2 COEFF')

%% Correlations of individual PCs and instant modulus

% corr(SCORE(:,1),y_train)
% corr(SCORE(:,2),y_train)
% corr(SCORE(:,3),y_train)
% corr(SCORE(:,4),y_train)
% corr(SCORE(:,5),y_train)

% not very good correlations individually.

%% for Principal Component Regression, we use more than one component.
% the example is copied from here:
% http://se.mathworks.com/help/stats/examples/partial-least-squares-regression-and-principal-components-regression.html
X = x_train;
y = y_train;
[n,p] = size(x_train);
[n_test,p_test] = size(x_test);

ncompmax = 50;

[PCALoadings,PCAScores,PCAVar] = princomp(X);
% in PCR, we combine the score values of multiple components.
% for example, the components 1-6
betaPCR = regress(y-mean(y), PCAScores(:,1:50));

betaPCR = PCALoadings(:,1:50)*betaPCR;
betaPCR = [mean(y) - mean(X)*betaPCR; betaPCR];

yfitPCR_train = [ones(n,1) x_train]*betaPCR;
yfitPCR_test = [ones(n_test,1) x_test]*betaPCR;

% a better correlation is achieved than using single components
corr(yfitPCR_train,y_train);
corr(yfitPCR_test,y_test);

r2train=0;
r2test=0;
rmsecvC=0;
rmsepC=0;

ncompmax = 15;

%4 key components to assess the performance of the model
r2train = zeros(ncompmax,1);
r2test = zeros(ncompmax,1);
rmsecvC = zeros(ncompmax,1);
rmsepC = zeros(ncompmax,1);
Error_Percent = zeros(ncompmax,1);

%Compute the 4 values for iterative addition of PLS components
h = waitbar(0);
for ncomp=1:ncompmax
    
    betaPCR = regress(y_train-mean(y_train), PCAScores(:,1:ncomp));
    betaPCR = PCALoadings(:,1:ncomp)*betaPCR;
    betaPCR = [mean(y_train) - mean(x_train)*betaPCR; betaPCR];
    
    yfit_PCR_train = [ones(n,1) x_train]*betaPCR;
    yfit_PCR_test = [ones(n_test,1) x_test]*betaPCR;
        
    s1 = sum((y_train - yfit_PCR_train).^2);
    s2 = sum((y_train - mean(y_train)).^2);
    r2train(ncomp) = (1 - (s1/s2))*100;
    rmsepC(ncomp) = sqrt(mean((y_test-yfit_PCR_test).^2));
    rmsecvC(ncomp) = sqrt(mean((y_train-yfit_PCR_train).^2));
    
    s11 = sum((y_test - yfit_PCR_test).^2);
    s22 = sum((y_test - mean(y_test)).^2);
    r2test(ncomp) = (1 - (s11/s22))*100;
    waitbar(ncomp / ncompmax);
    Error_Percent(ncomp) = (mean(abs(yfit_PCR_test-y_test)))/(max(y_train)-min(y_train))*100;
end
close(h)
% Plotting all the 4 parameters for quick comparison

figure;
title('Baseline per manuscript');
subplot(2,3,1);
plot(1:ncompmax,r2train)
xlabel('No.of PCA Comp');
ylabel('R^2 Train');
hold on
plot(find(r2train == min(r2train)),min(r2train),'o')
plot(find(r2test == max(r2test)),r2train,'X')
hold off

subplot(2,3,2);
plot(1:ncompmax,r2test)
xlabel('No.of PCA Comp');
ylabel('R^2 Test');
hold on
plot(find(r2test == max(r2test)),max(r2test),'o')
hold off

subplot(2,3,3);
plot(1:ncompmax,rmsepC)
xlabel('No.of PCA Comp');
ylabel('RMSEP');
hold on
plot(find(rmsepC == min(rmsepC)),min(rmsepC),'o')
plot(find(r2test == max(r2test)),rmsepC,'X')
hold off

subplot(2,3,4);
plot(1:ncompmax,rmsecvC)
xlabel('No.of PCA Comp');
ylabel('RMSECV');
hold on
plot(find(rmsecvC == min(rmsecvC)),min(rmsecvC),'o')
plot(find(r2test == max(r2test)),rmsecvC,'X')
hold off
title_to_display = strcat('PCA Regression :',func_param);
suptitle(title_to_display)

subplot(2,3,5);
plot(1:ncompmax,Error_Percent)
xlabel('No.of PCA Comp');
ylabel('Error %');
hold on
plot(find(Error_Percent == min(Error_Percent)),min(Error_Percent),'o')
plot(find(r2test == max(r2test)),Error_Percent,'X')
hold off

% Printing the values to respective matfile
C_pca = find(r2test == max(r2test)) ;
C = C_pca;
R2Train_pca = r2train(C);
R2Test_pca =r2test(C);
RMSECV_pca =rmsecvC(C);
RMSEP_pca =rmsepC(C);
Error_Percent_pca =Error_Percent(C);

if exist('D:\NIR Gui Project\results_regSelection\pca_record.mat','file') == 2
delete('D:\NIR Gui Project\results_regSelection\pca_record.mat');
else
save('D:\NIR Gui Project\results_regSelection\pca_record.mat','C_pca','R2Train_pca','RMSECV_pca','R2Test_pca','RMSEP_pca','Error_Percent_pca');
end

end

