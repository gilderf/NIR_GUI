function [ wavelength_final  ] = GA_var_processor( x_train, y_train, x_test , y_test )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% Select the data, X = spectra, Y = predicted feature
combination = [x_train y_train] ;

%test1 = gaplsopt(combination,1) ;
%test2 = gaplsopt(combination,2) ;

[b,c,d] = gaplssp(combination,182);

plotone(combination,b,c,d);
plotmore(combination,b);
plot(d);


%% TEST begins
A = 6;
[a,b] = sort(abs(d),'descend');

KFold = 10;
Indices = crossvalind('Kfold', numel(y_train), KFold);

% here we have additional loop (jj) which adds variables one by one based
% on the reliability index values

% jj starts from A (number of PLS components), because that is the minimum
% number of variables we need when we use A PLS components.
for jj=A:size(x_train,2)
    for ii=1:KFold

        ind_predict = find(Indices == ii);
        ind_train = find(Indices ~= ii);

        [XL,yl,XS,YS,beta,PCTVAR] = plsregress(x_train(ind_train,b(1:jj)),y_train(ind_train),A);
        yfit_GA_var(ind_predict,jj) = [ones(numel(ind_predict),1) x_train(ind_predict,b(1:jj))]*beta;

    end
end

%% Plot the correlation coefficient as a function of the number of variables
% figure
% plot(1:size(X,2),corr(yfit_GA_var,Y),'k','linewidth',2)
% hold on
% this find the maximum correlation
[corr_coefficient, var_no] = max(corr(yfit_GA_var,y_train));
% % mark the best correlation
% plot(var_no,corr_coefficient,'o')
% hold off
% xlabel('Number of variables','fontsize',18,'fontweight','bold')
% ylabel('Correlation coefficient','fontsize',18,'fontweight','bold')
% set(gca,'fontsize',16,'fontweight','bold')
% xlim([0 size(X,2)])
% 
% %% we can now plot the best model because we know how many variables are needed
% 
% yfit = yfit_GA_var(:,var_no);
% 
% figure
% scatter(Y,yfit)
% title(num2str(corr(Y,yfit)))


ga_vars = b(1: var_no);
wavelength_final = ga_vars ;

%[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,mcuve_vars),y_train,A);
ncompmax = 10;
%4 key components to assess the performance of the model
r2train = zeros(ncompmax,1);
r2test = zeros(ncompmax,1);
rmsecvC = zeros(ncompmax,1);
rmsepC = zeros(ncompmax,1);
Error_Percent = zeros(ncompmax,1);
ncomp = 0;

%Compute the 4 values for iterative addition of PLS components
for ncomp=1:ncompmax

%[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,ncomp,'CV',10);
[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,ga_vars),y_train,ncomp,'CV',10);
yfit_PLS_train = [ones(size(x_train,1),1) x_train(:,ga_vars)]*betaPLS;
yfit_PLS_test = [ones(size(x_test,1),1) x_test(:,ga_vars)]*betaPLS;


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

% Plotting all the 4 parameters for quick comparison
%Baseline setting

figure;
title('Baseline per manuscript');
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
suptitle('GA-PLS')

subplot(2,3,5);
plot(1:ncompmax,Error_Percent)
xlabel('No.of PLS Comp');
ylabel('Error %');
hold on
plot(find(Error_Percent == min(Error_Percent)),min(Error_Percent),'o')
plot(find(r2test == max(r2test)),Error_Percent,'X')
hold off


% Printing the values to respective matfile
C_gapls = find(r2test == max(r2test)) ;
C = C_gapls;
R2Train_gapls = r2train(C);
R2Test_gapls =r2test(C);
RMSECV_gapls =rmsecvC(C);
RMSEP_gapls =rmsepC(C);
Error_Percent_gapls =Error_Percent(C);


%if there is a previous copy of the file delete it and make a new one
if exist('D:\NIR Gui Project\results_varSelection\gapls_record.mat','file') == 2
delete('D:\NIR Gui Project\results_varSelection\gapls_record.mat');
end
save('D:\NIR Gui Project\results_varSelection\gapls_record.mat','C_gapls','R2Train_gapls','RMSECV_gapls','R2Test_gapls','RMSEP_gapls','Error_Percent_gapls');

end

