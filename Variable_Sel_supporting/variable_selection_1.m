%% Load and preprocess the data
% This files loads spectra and reference data, and does preprocessing to
% spectra.
% Needed files: data_average_spectrums.mat, Reference_data.xlsx
% NOTE: this is just copied from Jaakko's script. We may study also the
% preprocessing steps at some point.

% spectra are loaded, and wavelengths and spectra are stored in separate
% variables
myData = importdata('Data_average_spectrums.mat');
handles.wavelength = myData.Wavelength(895:1521);
handles.spectra = myData.Average_spectrums(:,895:1521);

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

%% Select the data, X = spectra, Y = predicted feature

X = Spectra_der2;
Y = handles.instant;

x_train =[X(1:260,:) ; X(281:528,:) ; X(554:653,:) ; X(679:869,:)];
x_test = [X(261:280,:) ;X(529:553,:) ;X(654:678,:)];

y_train =[Y(1:260,:) ; Y(281:528,:) ; Y(554:653,:) ; Y(679:869,:)];
y_test = [Y(261:280,:) ;Y(529:553,:) ;Y(654:678,:)];


%% MCUVE = Monte Carlo Uninformative Variable Elimination (PLS)
% In UVE, we estimate the "stability" of variables in the spectra
% for the PLS model. See better description from the articles.
% The principle is that we seek for variables that have a high stability,
% and use only those in our PLS model. 

% Number of PLS components (we do not know what is the optimal number,
% should be also studied)
A=3;
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

%% MCUVE
% How should we now select the variables based on the realiability index?
% One way is just to select a cutoff value (in this example 3) and select
% variables that have a higher index than 3.

cutoff_value = 3;
mcuve_vars = find(abs(UVE.RI) > cutoff_value);

% % then lets do a 10-fold cross-validation (we divide the data into 10
% % groups) and create PLS model using only variables indicated by mcuve_vars
% KFold = 10;
% Indices = crossvalind('Kfold', numel(y_train), KFold);
% 
% yfit_PLS_mcuve = zeros(size(y_train));
% for ii=1:KFold
%     % these values are used for training a sub-model
%     ind_train = find(Indices ~= ii);
%     % these values are predicted
%     ind_predict = find(Indices == ii);
% 
%     [XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(X(ind_train,mcuve_vars),Y(ind_train),A);
%     % predicted values are stored
%     yfit_PLS_mcuve(ind_predict) = [ones(numel(ind_predict),1) X(ind_predict,mcuve_vars)]*betaPLS;
% 
% end

[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,mcuve_vars),y_train,A);
%yfit_PLS_mcuve(ind_predict) = [ones(numel(ind_predict),1) X(ind_predict,mcuve_vars)]*betaPLS;


yfit_PLS_train_mcuve = [ones(size(x_train,1),1) x_train(:,mcuve_vars)]*betaPLS;
yfit_PLS_test_mcuve = [ones(size(x_test,1),1) x_test(:,mcuve_vars)]*betaPLS;


s1 = sum((y_train - yfit_PLS_train_mcuve).^2);
s2 = sum((y_train - mean(yfit_PLS_train_mcuve)).^2);
r2 = (1 - (s1/s2))*100
RMSEP = sqrt(mean((stats.Yresiduals).^2))
RMSECV = sqrt(msep(2,A))

s1 = sum((y_test - yfit_PLS_test_mcuve).^2);
s2 = sum((y_test - mean(yfit_PLS_test_mcuve)).^2);
r2 = (1 - (s1/s2))*100


% we can plot the result
figure
scatter(y_train,yfit_PLS_train_mcuve)
title([' trainset r = ' num2str(corr(y_train,yfit_PLS_train_mcuve))])

figure
scatter(y_test,yfit_PLS_test_mcuve)
title([' testset r = ' num2str(corr(y_test,yfit_PLS_test_mcuve))])

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
for jj=A:size(X,2)
    for ii=1:KFold

        ind_predict = find(Indices == ii);
        ind_train = find(Indices ~= ii);

        [XL,yl,XS,YS,beta,PCTVAR] = plsregress(x_train(ind_train,b(1:jj)),y_train(ind_train),A);
        yfit_PLS_uve_train(ind_predict,jj) = [ones(numel(ind_predict),1) x_train(ind_predict,b(1:jj))]*beta;
        %yfit_PLS_uve_test(ind_predict,jj) = [ones(numel(ind_predict),1) x_test(ind_predict,b(1:jj))]*beta;

    end
end

Indices = crossvalind('Kfold', numel(y_test), KFold);

% here we have additional loop (jj) which adds variables one by one based
% on the reliability index values

% jj starts from A (number of PLS components), because that is the minimum
% number of variables we need when we use A PLS components.
for jj=A:size(X,2)
    for ii=1:KFold

        ind_predict = find(Indices == ii);
        ind_train = find(Indices ~= ii);

        [XL,yl,XS,YS,beta,PCTVAR] = plsregress(x_train(ind_train,b(1:jj)),y_train(ind_train),A);
        %yfit_PLS_uve_train(ind_predict,jj) = [ones(numel(ind_predict),1) x_train(ind_predict,b(1:jj))]*beta;
        yfit_PLS_uve_test(ind_predict,jj) = [ones(numel(ind_predict),1) x_test(ind_predict,b(1:jj))]*beta;

    end
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

%% Plot the correlation coefficient as a function of the number of variables
figure
plot(1:size(x_test,2),corr(yfit_PLS_uve_test,y_test),'k','linewidth',2)
hold on
% this find the maximum correlation
[corr_coefficient, var_no2] = max(corr(yfit_PLS_uve_test,y_test));
% mark the best correlation
plot(var_no2,corr_coefficient,'o')
hold off
xlabel('Number of variables','fontsize',18,'fontweight','bold')
ylabel('Correlation coefficient','fontsize',18,'fontweight','bold')
set(gca,'fontsize',16,'fontweight','bold')
xlim([0 size(x_test,2)])

%% we can now plot the best model because we know how many variables are needed

yfit = yfit_PLS_uve_train(:,var_no);

figure
scatter(y_train,yfit)
title(num2str(corr(y_train,yfit)))

yfit2 = yfit_PLS_uve_test(:,var_no2);

figure
scatter(y_test,yfit2)
title(num2str(corr(y_test,yfit2)))


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


% lets create a final model based on the selected variables

KFold = 10;
Indices = crossvalind('Kfold', numel(y_train), KFold);
yfit_PLS_CARS_train = zeros(size(y_train));

% Cross-validation is conducted similarly as earlier
for ii=1:KFold
    ind_predict = find(Indices == ii);
    ind_train = find(Indices ~= ii);

    [XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(ind_train,CARS.vsel),y_train(ind_train),CARS.optLV);
    yfit_PLS_CARS_train(ind_predict) = [ones(numel(ind_predict),1) x_train(ind_predict,CARS.vsel)]*betaPLS;
end

%% Here we can plot the final model of CARS
figure
scatter(y_train,yfit_PLS_CARS_train)
%title(num2str(corr(y_train,yfit_PLS_CARS_train)))
title([' trainset r = ' num2str(corr(y_train,yfit_PLS_CARS_train))])

%% TEST DATA

% lets create a final model based on the selected variables

KFold = 10;
Indices = crossvalind('Kfold', numel(y_test), KFold);
yfit_PLS_CARS_test = zeros(size(y_test));

% Cross-validation is conducted similarly as earlier
for ii=1:KFold
    ind_predict = find(Indices == ii);
    ind_train = find(Indices ~= ii);

    [XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_test(ind_train,CARS.vsel),y_test(ind_train),CARS.optLV);
    yfit_PLS_CARS_test(ind_predict) = [ones(numel(ind_predict),1) x_test(ind_predict,CARS.vsel)]*betaPLS;
end

%% Here we can plot the final model of CARS
figure
scatter(y_test,yfit_PLS_CARS_test)
%title(num2str(corr(y_test,yfit_PLS_CARS_test)))
title([' testset r = ' num2str(corr(y_test,yfit_PLS_CARS_test))])

s1 = sum((y_train - yfit_PLS_CARS_train).^2);
s2 = sum((y_train - mean(yfit_PLS_CARS_train)).^2);
r2 = (1 - (s1/s2))*100
RMSEP = sqrt(mean((stats.Yresiduals).^2))
RMSECV = sqrt(msep(2,A))

s1 = sum((y_test - yfit_PLS_CARS_test).^2);
s2 = sum((y_test - mean(yfit_PLS_CARS_test)).^2);
r2 = (1 - (s1/s2))*100

