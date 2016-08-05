function [ wavelength_final  ] = Jack_Knife_var_processor( x_train, y_train, x_test , y_test)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

testingJK=jackknife_vs(x_train,y_train,6);

selected_variables=testingJK;
wavelength_final = selected_variables ;

ncompmax = 15;
%4 key components to assess the performance of the model
r2train = zeros(ncompmax,1);
r2test = zeros(ncompmax,1);
rmsecvC = zeros(ncompmax,1);
rmsepC = zeros(ncompmax,1);
Error_Percent = zeros(ncompmax,1);
ncomp = 0;

h = waitbar(0);
%Compute the 4 values for iterative addition of PLS components
for ncomp=1:ncompmax

%[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,ncomp,'CV',10);
[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,selected_variables),y_train,ncomp,'CV',10);
yfit_PLS_train = [ones(size(x_train,1),1) x_train(:,selected_variables)]*betaPLS;
yfit_PLS_test = [ones(size(x_test,1),1) x_test(:,selected_variables)]*betaPLS;

waitbar(ncomp / ncompmax);

s1 = sum((y_train - yfit_PLS_train).^2);
s2 = sum((y_train - mean(y_train)).^2);
r2train(ncomp) = (1 - (s1/s2))*100;
rmsepC(ncomp) = sqrt(mean((y_test - yfit_PLS_test).^2));
rmsecvC(ncomp) = sqrt(msep(2,ncomp+1));

s11 = sum((y_test - yfit_PLS_test).^2);
s22 = sum((y_test - mean(y_test)).^2);
r2test(ncomp) = (1 - (s11/s22))*100;

Error_Percent(ncomp) = abs(mean(abs(yfit_PLS_test-y_test)))/(max(y_train)- min(y_train))*100;
end
close(h)
% Plotting all the 4 parameters for quick comparison
%Baseline setting

figure;
%title('Baseline per manuscript');
subplot(2,3,1);
plot(1:ncompmax,r2train)
xlabel('No.of PLS Comp');
ylabel('R^2 Train');
hold on
plot(find(r2train == min(r2train)),min(r2train),'o')
plot(find(r2test == max(r2test)),r2train,'X')
hold off

subplot(2,3,3);
plot(1:ncompmax,r2test)
xlabel('No.of PLS Comp');
ylabel('R^2 Test');
hold on
plot(find(r2test == max(r2test)),max(r2test),'o')
hold off

subplot(2,3,4);
plot(1:ncompmax,rmsepC)
xlabel('No.of PLS Comp');
ylabel('RMSEP');
hold on
plot(find(rmsepC == min(rmsepC)),min(rmsepC),'o')
plot(find(r2test == max(r2test)),rmsepC,'X')
hold off

subplot(2,3,2);
plot(1:ncompmax,rmsecvC)
xlabel('No.of PLS Comp');
ylabel('RMSECV');
hold on
plot(find(rmsecvC == min(rmsecvC)),min(rmsecvC),'o')
plot(find(r2test == max(r2test)),rmsecvC,'X')
hold off
suptitle('Jackknife-PLS')

subplot(2,3,5);
plot(1:ncompmax,Error_Percent)
xlabel('No.of PLS Comp');
ylabel('Error %');
hold on
plot(find(Error_Percent == min(Error_Percent)),min(Error_Percent),'o')
plot(find(r2test == max(r2test)),Error_Percent,'X')
hold off


% Printing the values to respective matfile
C_jkpls = find(r2test == max(r2test)) ;
C = C_jkpls;
R2Train_jkpls = r2train(C);
R2Test_jkpls =r2test(C);
RMSECV_jkpls =rmsecvC(C);
RMSEP_jkpls =rmsepC(C);
Error_Percent_jkpls =Error_Percent(C);

%if there is a previous copy of the file delete it and make a new one
if exist('D:\NIR Gui Project\results_varSelection\jkpls_record.mat','file') == 2
delete('D:\NIR Gui Project\results_varSelection\jkpls_record.mat');
end
save('D:\NIR Gui Project\results_varSelection\jkpls_record.mat','C_jkpls','R2Train_jkpls','RMSECV_jkpls','R2Test_jkpls','RMSEP_jkpls','Error_Percent_jkpls');

end

