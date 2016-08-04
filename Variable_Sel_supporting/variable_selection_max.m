%% Load and preprocess the data
ccc;
myData = importdata('Data_average_spectrums.mat');
handles.wavelength = myData.Wavelength(1:2000);
handles.spectra = myData.Average_spectrums(:,1:2000);

% Reference data from the excel file are loaded and stored in separate
% variables
myRef = xlsread('Reference_data.xlsx');
handles.thickness = myRef(:,1);
handles.thickness_cal = myRef(:,2);
handles.ICRS = myRef(:,3);
handles.instant = myRef(:,5);
handles.equilibrium = myRef(:,7);
handles.freq1 = myRef(:,8);
handles.freq2 = myRef(:,9);
handles.freq3 = myRef(:,10);
handles.freq4 = myRef(:,11);
handles.freq5 = myRef(:,12);
handles.freq6 = myRef(:,13);
handles.freq7 = myRef(:,14);
clear fn path myRef myData

% This does preprocessing to spectra. 
% First, Savitzky-Golay smoothing is applied. If you are not yet familiar 
% with Savitzky-Golay, you should read a little bit about it. It is quite 
% simple, and a very common method in spectroscopy.

% Second, the second derivative spectra area calculated

% Third, a moving-average filter is applied to reduce noise.

% Best smoothing    1) window 40, degree 3
%                   2) Second derivative
%                   3) window 10, degree 5
for kk = 1:size(handles.spectra,1)
    Spectra_smooth = smooth(handles.wavelength,handles.spectra(kk,:),40,'sgolay',3);
    Spectra_smooth = Spectra_smooth';
    for ll = 1:size(Spectra_smooth,2)-2
        Spectra_der(:,ll) = (Spectra_smooth(:,ll+2)-2*Spectra_smooth(:,ll+1)+Spectra_smooth(:,ll))./(((handles.wavelength(ll+2)-handles.wavelength(ll))/2).^2);
    end
    wind = 10; % Normal 10
    Spectra_der2(kk,:) = filter(ones(1,wind)/wind,1,Spectra_der);
end
clear Spectra_smooth kk ll  der Spectra_der

% % The best test set was found to be Sample 10B, Sample 10I and Sample 5G
% ind_1 = 529; % Sample 10B - ICRS 0
% ind_2 = 553;
% ind_3 = 654; % Sample 10I - ICRS 1
% ind_4 = 678;
% ind_5 = 261; % Sample 5G - ICRS 2
% ind_6 = 280;

% Wavelength indexes
wave_start = 940;
wave_end = 1340;

%% Lets use variable X for spectra and variable Y for predicted data.
% In this case, lets use second derivative spectra (we could use the
% original spectra without preprocessing as well)
X = Spectra_der2(:,wave_start:wave_end);
x_train =[X(1:260,:) ; X(281:528,:) ; X(554:653,:) ; X(679:869,:)];
x_test = [X(261:280,:) ;X(529:553,:) ;X(654:678,:)];

% And the predicted feature is instant modulus
% Y = handles.instant;
% Y = Y./1e6;
Y = handles.equilibrium;
y_train =[Y(1:260,:) ; Y(281:528,:) ; Y(554:653,:) ; Y(679:869,:)];
y_test = [Y(261:280,:) ;Y(529:553,:) ;Y(654:678,:)];

x_train = x_train(find(~isnan(y_train)),:);
y_train = y_train(find(~isnan(y_train)));

x_test = x_test(find(~isnan(y_test)),:);
y_test = y_test(find(~isnan(y_test)));


%% MCUVE = Monte Carlo Uninformative Variable Elimination (PLS)
% In UVE, we estimate the "stability" of variables in the spectra
% for the PLS model. See better description from the articles.
% The principle is that we seek for variables that have a high stability,
% and use only those in our PLS model. 

% Number of PLS components (we do not know what is the optimal number,
% should be also studied)
A=6;
% Pre-processing method (I think the function requires some preprocessing
% method to be used. 'center' just centers the mean of the spectrum to 0,
% which in case of 2nd derivative spectrum is already the case (at least
% close to it). 
method='center';
% Parameters related to MC estimation.
N=1000;
ratio=0.75;
UVE=mcuvepls(x_train,y_train,A,method,N,ratio);

% UVE.RI shows us the reliability of the variables. 
figure
plot(abs(UVE.RI),'linewidth',2);
xlabel('variable index');
ylabel('reliability index');
set(gcf,'color','w');

%% MCUVE, continued
% it is difficult to know what is the best cutoff value. We can also just
% add variables one by one to the model based on reliability index and
% observe how the model prediction progesses.

% lets sort the RI values in descending order
[a,b] = sort(abs(UVE.RI),'descend');

KFold = 10;
Indices = crossvalind('Kfold', numel(y_train), KFold);

% here we have additional loop (jj) which adds variables one by one based
% on the reliability index values

% jj starts from A (number of PLS components), because that is the minimum
% number of variables we need when we use A PLS components.
for jj=A:size(x_train,2)
    
        [XL,yl,XS,YS,beta,PCTVAR] = plsregress(x_train(:,b(1:jj)),y_train,A);
        yfit_PLS_uve_train(:,jj) = [ones(size(x_train,1),1) x_train(:,b(1:jj))]*beta;        
end

%% Plot the correlation coefficient as a function of the number of variables
figure
plot(1:size(x_train,2),corr(yfit_PLS_uve_train,y_train),'k','linewidth',2)
hold on
% this find the maximum correlation
[corr_coefficient, var_no] = max(corr(yfit_PLS_uve_train,y_train));
% mark the best correlation
plot(var_no,corr_coefficient,'o')
hold off
xlabel('Number of variables','fontsize',18,'fontweight','bold')
ylabel('Correlation coefficient','fontsize',18,'fontweight','bold')
set(gca,'fontsize',16,'fontweight','bold')
xlim([0 size(x_train,2)])

%% we can now plot the best model because we know how many variables are needed

yfit_train = yfit_PLS_uve_train(:,var_no);

figure
scatter(y_train,yfit_train)
title(num2str(corr(y_train,yfit_train)))

mcuve_vars = b(1: var_no);

%[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,mcuve_vars),y_train,A);
ncompmax = 15;
%4 key components to assess the performance of the model
r2train = zeros(ncompmax,1);
r2test = zeros(ncompmax,1);
rmsecvC = zeros(ncompmax,1);
rmsepC = zeros(ncompmax,1);
ncomp = 0;

%Compute the 4 values for iterative addition of PLS components
for ncomp=1:ncompmax

%[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,ncomp,'CV',10);
[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,mcuve_vars),y_train,ncomp,'CV',10);
yfit_PLS_train = [ones(size(x_train,1),1) x_train(:,mcuve_vars)]*betaPLS;
yfit_PLS_test = [ones(size(x_test,1),1) x_test(:,mcuve_vars)]*betaPLS;


s1 = sum((y_train - yfit_PLS_train).^2);
s2 = sum((y_train - mean(y_train)).^2);
r2train(ncomp) = (1 - (s1/s2))*100;
rmsepC(ncomp) = sqrt(mean((stats.Yresiduals).^2));
rmsecvC(ncomp) = sqrt(msep(2,ncomp+1));

s11 = sum((y_test - yfit_PLS_test).^2);
s22 = sum((y_test - mean(y_test)).^2);
r2test(ncomp) = (1 - (s11/s22))*100;
end

% Plotting all the 4 parameters for quick comparison
%Baseline setting

figure;
title('Baseline per manuscript');
subplot(2,2,1);
plot(1:ncompmax,r2train)
xlabel('No.of PLS Comp');
ylabel('R^2 Train');
hold on
plot(find(r2train == min(r2train)),min(r2train),'o')
plot(find(r2test == max(r2test)),r2train,'X')
hold off

subplot(2,2,3);
plot(1:ncompmax,r2test)
xlabel('No.of PLS Comp');
ylabel('R^2 Test');
hold on
plot(find(r2test == max(r2test)),max(r2test),'o')
hold off

subplot(2,2,4);
plot(1:ncompmax,rmsepC)
xlabel('No.of PLS Comp');
ylabel('RMSEP');
hold on
plot(find(rmsepC == min(rmsepC)),min(rmsepC),'o')
plot(find(r2test == max(r2test)),rmsepC,'X')
hold off

subplot(2,2,2);
plot(1:ncompmax,rmsecvC)
xlabel('No.of PLS Comp');
ylabel('RMSECV');
hold on
plot(find(rmsecvC == min(rmsecvC)),min(rmsecvC),'o')
plot(find(r2test == max(r2test)),rmsecvC,'X')
hold off
suptitle('UVE-PLS')


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

CARS.vsel;

%[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,mcuve_vars),y_train,A);
ncompmax = 15;
%4 key components to assess the performance of the model
r2train = zeros(ncompmax,1);
r2test = zeros(ncompmax,1);
rmsecvC = zeros(ncompmax,1);
rmsepC = zeros(ncompmax,1);
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
rmsepC(ncomp) = sqrt(mean((stats.Yresiduals).^2));
rmsecvC(ncomp) = sqrt(msep(2,ncomp+1));

s11 = sum((y_test - yfit_PLS_test).^2);
s22 = sum((y_test - mean(y_test)).^2);
r2test(ncomp) = (1 - (s11/s22))*100;
end

% Plotting all the 4 parameters for quick comparison
%Baseline setting

figure;
%title('Baseline per manuscript');
subplot(2,2,1);
plot(1:ncompmax,r2train)
xlabel('No.of PLS Comp');
ylabel('R^2 Train');
hold on
plot(find(r2train == min(r2train)),min(r2train),'o')
plot(find(r2test == max(r2test)),r2train,'X')
hold off

subplot(2,2,3);
plot(1:ncompmax,r2test)
xlabel('No.of PLS Comp');
ylabel('R^2 Test');
hold on
plot(find(r2test == max(r2test)),max(r2test),'o')
hold off

subplot(2,2,4);
plot(1:ncompmax,rmsepC)
xlabel('No.of PLS Comp');
ylabel('RMSEP');
hold on
plot(find(rmsepC == min(rmsepC)),min(rmsepC),'o')
plot(find(r2test == max(r2test)),rmsepC,'X')
hold off

subplot(2,2,2);
plot(1:ncompmax,rmsecvC)
xlabel('No.of PLS Comp');
ylabel('RMSECV');
hold on
plot(find(rmsecvC == min(rmsecvC)),min(rmsecvC),'o')
plot(find(r2test == max(r2test)),rmsecvC,'X')
hold off
suptitle('CARS-PLS')



