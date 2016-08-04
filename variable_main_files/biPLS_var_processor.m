function [ wavelength_final  ] = biPLS_var_processor( x_train, y_train, x_test , y_test ,range )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


xaxis = range;
% Y = handles.instant;

biModel=bipls(x_train,y_train,10,'mean',20,xaxis,'syst123',5);
biplstable(biModel); 
selected_intervals = [5 17 19] ;
FinalModel=plsmodel(biModel,selected_intervals,10,'mean','syst123',3); 
% plsrmse(FinalModel);
% plspvsm(FinalModel,6);
% 
% predModel=plspredict(x_test,FinalModel,6,y_test);
% plsrmse(predModel);
% plspvsm(predModel,6);

bipls_var = [];
for i=1:max(size(FinalModel.selected_intervals))
        temp=FinalModel.allint(FinalModel.selected_intervals(i),2):FinalModel.allint(FinalModel.selected_intervals(i),3);
        bipls_var=[bipls_var temp];
    
end

wavelength_final = bipls_var;

%[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,mcuve_vars),y_train,A);
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

%[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,ncomp,'CV',10);
[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,bipls_var),y_train,ncomp,'CV',10);
yfit_PLS_train = [ones(size(x_train,1),1) x_train(:,bipls_var)]*betaPLS;
yfit_PLS_test = [ones(size(x_test,1),1) x_test(:,bipls_var)]*betaPLS;


s1 = sum((y_train - yfit_PLS_train).^2);
s2 = sum((y_train - mean(y_train)).^2);
r2train(ncomp) = (1 - (s1/s2))*100;
rmsepC(ncomp) = sqrt(mean((y_test - yfit_PLS_test).^2));
rmsecvC(ncomp) = sqrt(mean((stats.Yresiduals).^2));

s11 = sum((y_test - yfit_PLS_test).^2);
s22 = sum((y_test - mean(y_test)).^2);
r2test(ncomp) = (1 - (s11/s22))*100;
waitbar(ncomp / ncompmax);

Error_Percent(ncomp) = (mean(abs(yfit_PLS_test-y_test)))/(max(y_train)- min(y_train))*100;

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
suptitle('biPLS')

subplot(2,3,5);
plot(1:ncompmax,Error_Percent)
xlabel('No.of PLS Comp');
ylabel('Error %');
hold on
plot(find(Error_Percent == min(Error_Percent)),min(Error_Percent),'o')
plot(find(r2test == max(r2test)),Error_Percent,'X')
hold off

% Printing the values to respective matfile
C_bipls = find(r2test == max(r2test)) ;
C = C_bipls;
R2Train_bipls = r2train(C);
R2Test_bipls =r2test(C);
RMSECV_bipls =rmsecvC(C);
RMSEP_bipls =rmsepC(C);
Error_Percent_bipls =Error_Percent(C);

if exist('D:\NIR Gui Project\results_varSelection\bipls_record.mat','file') == 2
delete('D:\NIR Gui Project\results_varSelection\bipls_record.mat');
else
save('D:\NIR Gui Project\results_varSelection\bipls_record.mat','C_bipls','R2Train_bipls','RMSECV_bipls','R2Test_bipls','RMSEP_bipls','Error_Percent_bipls');
end

end

