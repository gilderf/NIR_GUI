function [ r2train,r2test,rmsecvC,rmsepC  ] = lasso_regression( x_train, y_train, x_test , y_test, func_param )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
tic
h = waitbar(0);
[B,S] = lasso(x_train,y_train,'CV',10);
close(h)
toc
Xplusn_test = [ones(size(x_test,1),1) x_test];
Xplusn_train = [ones(size(x_train,1),1) x_train];

limitSD = size(S.DF,2);
r2train = zeros(limitSD,1);
rmsecvC = zeros(limitSD,1);
r2test = zeros(limitSD,1);
rmsepC = zeros(limitSD,1);
Error_Percent = zeros(limitSD,1);

for n = 1 : size(S.DF,2)
    
    fitDF61n_train = Xplusn_train * [S.Intercept(n); B(:,n)];
    fitDF61n_test = Xplusn_test * [S.Intercept(n); B(:,n)];
    
    s1 = sum((y_train - fitDF61n_train).^2);
    s2 = sum((y_train - mean(y_train)).^2);
    r2train(n) = (1 - (s1/s2))*100;
    rmsecvC(n) = sqrt(mean((y_train-fitDF61n_train).^2));
    
    s11 = sum((y_test - fitDF61n_test).^2);
    s22 = sum((y_test - mean(y_test)).^2);
    r2test(n) = (1 - (s11/s22))*100;
    rmsepC(n) = sqrt(mean((y_test - fitDF61n_test).^2));
    
    Error_Percent(n) = abs(mean(abs(fitDF61n_test-y_test)))/(max(y_train)-min(y_train))*100;

end

% Print the results Grid
ncompmax = size(S.DF,2);
figure;
title('Baseline per manuscript');
subplot(2,3,1);
plot(1:ncompmax,r2train)
xlabel('Lasso Index');
ylabel('R^2 Train');
hold on
plot(find(r2train == min(r2train)),min(r2train),'o')
plot(find(r2test == max(r2test)),r2train,'X')
hold off

subplot(2,3,3);
plot(1:ncompmax,r2test)
xlabel('Lasso Index');
ylabel('R^2 Test');
hold on
plot(find(r2test == max(r2test)),max(r2test),'o')
hold off

subplot(2,3,4);
plot(1:ncompmax,rmsepC)
xlabel('Lasso Index');
ylabel('rmsepC');
hold on
plot(find(rmsepC == min(rmsepC)),min(rmsepC),'o')
plot(find(r2test == max(r2test)),rmsepC,'X')
hold off

subplot(2,3,2);
plot(1:ncompmax,rmsecvC)
xlabel('Lasso Index');
ylabel('RMSECV');
hold on
plot(find(rmsecvC == min(rmsecvC)),min(rmsecvC),'o')
plot(find(r2test == max(r2test)),rmsecvC,'X')
hold off

title_to_display = strcat('LASSO Regression : ',func_param);
suptitle(title_to_display)

subplot(2,3,5);
plot(1:ncompmax,Error_Percent)
xlabel('Lambda Value');
ylabel('Error %');
hold on
plot(find(Error_Percent == min(Error_Percent)),min(Error_Percent),'o')
plot(find(r2test == max(r2test)),Error_Percent,'X')
hold off

% Printing the values to respective matfile
C_lasso = find(r2test == max(r2test)) ;
C = C_lasso;
R2Train_lasso = r2train(C);
R2Test_lasso = r2test(C);
RMSECV_lasso = rmsecvC(C);
RMSEP_lasso = rmsepC(C);
Error_Percent_lasso = Error_Percent(C);

%if there is a previous copy of the file delete it and make a new one
if exist('D:\NIR Gui Project\results_regSelection\lasso_record.mat','file') == 2
delete('D:\NIR Gui Project\results_regSelection\lasso_record.mat');
end
save('D:\NIR Gui Project\results_regSelection\lasso_record.mat','C_lasso','R2Train_lasso','RMSECV_lasso','R2Test_lasso','RMSEP_lasso','Error_Percent_lasso');

end

