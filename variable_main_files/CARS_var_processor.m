function [ wavelength_final  ] = CARS_var_processor( x_train, y_train, x_test , y_test )
%% Competitive adaptive reweighted sampling (CARS)
% lets move on to the next method, CARS. Again, see the articles (and my
% manuscript) for better description.

% X = spectra 
% Y = predicted feature 
% 10 = maximum number of PLS components (function automatically find the
% best between 1-10 if we set 10 to be the maximum)
% 10 = K-Fold in cross-validation
% 'center' = pre-processing method
% 100 = number of sampling runs in Monte Carlo
% 0 = Number of PLS components is chosen based on global minimum of
% prediction error
% 1 = 'original version' = the regression coefficients are used as the likelihood
% values for the coefficients to be used (they are not just sorted in
% descending order and selected based on that. If we want that to be used,
% then we select = 0) 
% see original article.

CARS=carspls(x_train,y_train,10,10,'center',100,0,1);

% CARS contains some useful information such as CARS.optLV = number of PLS
% components, and CARS.vsel = selected variables

% This plots one type of presenation of the results, the original article
% explains what is going on here...
figure, plotcars(CARS)

wavelength_final = CARS.vsel;

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
for ncomp=1:ncompmax

%[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,ncomp,'CV',10);
[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,CARS.vsel),y_train,ncomp,'CV',10);
yfit_PLS_train = [ones(size(x_train,1),1) x_train(:,CARS.vsel)]*betaPLS;
yfit_PLS_test = [ones(size(x_test,1),1) x_test(:,CARS.vsel)]*betaPLS;


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
suptitle('CARS-PLS')

subplot(2,3,5);
plot(1:ncompmax,Error_Percent)
xlabel('No.of PLS Comp');
ylabel('Error %');
hold on
plot(find(Error_Percent == min(Error_Percent)),min(Error_Percent),'o')
plot(find(r2test == max(r2test)),Error_Percent,'X')
hold off

% Printing the values to respective matfile
C_cars = find(r2test == max(r2test)) ;
C = C_cars;
R2Train_cars = r2train(C);
R2Test_cars =r2test(C);
RMSECV_cars =rmsecvC(C);
RMSEP_cars =rmsepC(C);
Error_Percent_cars =Error_Percent(C);

%if there is a previous copy of the file delete it and make a new one
if exist('D:\NIR Gui Project\results_varSelection\cars_record.mat','file') == 2
delete('D:\NIR Gui Project\results_varSelection\cars_record.mat');
end
save('D:\NIR Gui Project\results_varSelection\cars_record.mat','C_cars','R2Train_cars','RMSECV_cars','R2Test_cars','RMSEP_cars','Error_Percent_cars');

end

