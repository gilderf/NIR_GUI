function [ r2train,r2test,rmsecvC,rmsepC  ] = pls_regression( x_train, y_train, x_test , y_test, func_param )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


%% optimum number of PLS components are still not known
%Lets compute the till 15 PLS components
ncompmax = 15;
%4 key components to assess the performance of the model
r2train = zeros(ncompmax,1);
r2test = zeros(ncompmax,1);
rmsecvC = zeros(ncompmax,1);
rmsepC = zeros(ncompmax,1);
Error_Percent = zeros(ncompmax,1);
ncomp = 0;

%Compute the 4 values for iterative addition of PLS components
h = waitbar(0);
for ncomp=1:ncompmax

[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,ncomp,'CV',10);
yfit_PLS_train = [ones(size(x_train,1),1) x_train]*betaPLS;
yfit_PLS_test = [ones(size(x_test,1),1) x_test]*betaPLS;


s1 = sum((y_train - yfit_PLS_train).^2);
s2 = sum((y_train - mean(y_train)).^2);
r2train(ncomp) = (1 - (s1/s2))*100;
rmsepC(ncomp) = sqrt(mean((y_test - yfit_PLS_test).^2));
rmsecvC(ncomp) = sqrt(mean((stats.Yresiduals).^2));

s11 = sum((y_test - yfit_PLS_test).^2);
s22 = sum((y_test - mean(y_test)).^2);
r2test(ncomp) = (1 - (s11/s22))*100;
waitbar(ncomp / ncompmax);
Error_Percent(ncomp) = (mean(abs(yfit_PLS_test-y_test)))/(max(y_train)-min(y_train))*100;
end
close(h)

% Plotting all the 4 parameters for quick comparison
%Baseline setting

figure;
title('PLS');
subplot(2,3,1);
plot(1:ncompmax,r2train);
xlabel('No.of PLS Comp');
ylabel('R^2 Train');
hold on
plot(find(r2train == min(r2train)),min(r2train),'o');
plot(find(r2test == max(r2test)),r2train,'X');
hold off

subplot(2,3,2);
plot(1:ncompmax,r2test);
xlabel('No.of PLS Comp');
ylabel('R^2 Test');
hold on
plot(find(r2test == max(r2test)),max(r2test),'o');
hold off

subplot(2,3,3);
plot(1:ncompmax,rmsepC);
xlabel('No.of PLS Comp');
ylabel('RMSEP');
hold on
plot(find(rmsepC == min(rmsepC)),min(rmsepC),'o');
plot(find(r2test == max(r2test)),rmsepC,'X');
hold off

subplot(2,3,4);
plot(1:ncompmax,rmsecvC);
xlabel('No.of PLS Comp');
ylabel('RMSECV');
hold on
plot(find(rmsecvC == min(rmsecvC)),min(rmsecvC),'o');
plot(find(r2test == max(r2test)),rmsecvC,'X');
hold off
title_to_display = strcat('PLS Regression :',func_param);
suptitle(title_to_display)

subplot(2,3,5);
plot(1:ncompmax,Error_Percent);
xlabel('No.of PLS Comp');
ylabel('Error %');
hold on
plot(find(Error_Percent == min(Error_Percent)),min(Error_Percent),'o');
plot(find(r2test == max(r2test)),Error_Percent,'X');
hold off

% Printing the values to respective matfile
C_pls = find(r2test == max(r2test)) ;
C = C_pls;
R2Train_pls = r2train(C);
R2Test_pls =r2test(C);
RMSECV_pls =rmsecvC(C);
RMSEP_pls =rmsepC(C);
Error_Percent_pls =Error_Percent(C);

%if there is a previous copy of the file delete it and make a new one
if exist('D:\NIR Gui Project\results_regSelection\pls_record.mat','file') == 2
delete('D:\NIR Gui Project\results_regSelection\pls_record.mat');
end
save('D:\NIR Gui Project\results_regSelection\pls_record.mat','C_pls','R2Train_pls','RMSECV_pls','R2Test_pls','RMSEP_pls','Error_Percent_pls');

end

