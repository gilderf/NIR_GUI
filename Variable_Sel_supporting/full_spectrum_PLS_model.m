%% Load and preprocess the data
% This files loads spectra and reference data, and does preprocessing to
% spectra.
% Needed files: data_average_spectrums.mat, Reference_data.xlsx
ccc

% spectra are loaded, and wavelengths and spectra are stored in separate
% variables
myData = importdata('Data_average_spectrums.mat');
handles.wavelength = myData.Wavelength(1:2000);
handles.spectra = myData.Average_spectrums(:,1:2000);
% handles.wavelength = myData.Wavelength(895:1521);
% handles.spectra = myData.Average_spectrums(:,895:1521);
% handles.wavelength = myData.Wavelength(940:1340);
% handles.spectra = myData.Average_spectrums(:,940:1340);

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
% The best test set was found to be Sample 10B, Sample 10I and Sample 5G
ind_1 = 529; % Sample 10B - ICRS 0
ind_2 = 553;
ind_3 = 654; % Sample 10I - ICRS 1
ind_4 = 678;
ind_5 = 261; % Sample 5G - ICRS 2
ind_6 = 280;

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
Y = handles.instant;
Y = Y./1e6;
y_train =[Y(1:260,:) ; Y(281:528,:) ; Y(554:653,:) ; Y(679:869,:)];
y_test = [Y(261:280,:) ;Y(529:553,:) ;Y(654:678,:)];

%% PLS regression, basics
% So we know some basics already. PLSREGRESS function is easy to use, we
% just give it the spectra (X), the predicted feature (Y) and number of PLS
% components (e.g. 10). Function returns us the "PLS Component spectra",
% the regression coefficients for each variable (betaPLS), percentage of
% explained covariance by each variable (PctVar), prediction error (msep)
% and some additional information in stats.
[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,5);

% now we could use the model to predict instant modulus from the spectra by
% using betaPLS:
yfit = [ones(size(x_test,1),1) x_test]*betaPLS;

% we can plot a scatter plot and calculate the correlation between the
% reference values and predicted values
figure
scatter(y_test,yfit)
title(['r = ' num2str(corr(y_test,yfit))])

%%
% there is one problem: we are predicting the values of the samples that
% were already used when building the model. We should have an external
% test sample set, or then we can use cross-validation (if you are not
% familiar with it, read about it). In cross-validation, a portion of
% samples acts as validation samples, and this is repeated so that every
% sample acts as validation samples. This means that during
% cross-validation, multiple PLS models are built, and the results of
% predictions of each run are stored and in the combined to evaluate the
% performance of the final model.

% First, we can use PLSREGRESS-functions own cross-validation when we
% study, how many PLS components we should use in our model.

ncomp_max = 15;

% 'cv',10 option tells the function to conduct 10-fold
% cross-validation, which means that the data is divided into 10
% groups. So, 10 pls regression models are created, and from each run
% the predictions of the validation group are investigated. Based on
% these, we can calculate the cross-validated prediction error
% (root-mean-square error of cross validation, or RMSECV). We do this
% for different number of PLS components, in this example, 1-15. When
% we plot the RMSECV as a function of the number of PLS components, we
% can see what would be the optimal number of PLS components.

rmsecv = zeros(ncomp_max,1);
for ncomp=1:ncomp_max
    [XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,ncomp,'cv',10);
    rmsecv(ncomp) = sqrt(msep(2,ncomp+1));
end

figure
plot(1:ncomp_max,rmsecv)
hold on
plot(find(rmsecv == min(rmsecv)),min(rmsecv),'o')
hold off

% from the figure, we can see that RMSECV first decreases as we increase
% the number of PLS components. However, after 6 PLS components, RMSECV
% actually starts to increase! This is an "ideal" RMSECV curve in that
% sense that it is very easy to select the optimal number of PLS components
% based on the minimum of RMSECV: we should use 6 PLS components. Sometimes
% this kind of clear minimum is not observed, and selection gets more
% difficult. Always a simpler model is preferred, so if there is only very
% small decrease in RMSECV as the components are added, we should select
% the simpler model (=less PLS components).
%
% It should be noted that this particular RMSECV curve applies only to
% these models that are based on "full spectrum" (=we haven't done any kind
% of variable selection) second derivative spectra. If we do some variable
% selection or even if we use the "raw" spectra instead of second
% derivative spectra, the curve may become different.

%% Now that we know the optimal number of PLS components, we can build
% an cross-validated PLS model.

ncomp = 6;
% we first initiate a vector that is used to store the predictions from
% cross-validation. The size is obviously equal to the size of predicted
% samples.
%yfit_PLS = zeros(size(y_test));

% For cross-validation, we need to select the "Fold" of cross-validation,
% which means to how many groups we divide our samples in cross-validation.
% KFold = 10;
% 
% % crossvalind creates indices which indicate to which group each sample
% % belongs. The sizes of these groups are equal if possible (depends on if
% % the number of samples can be equally divided into chosen number of sample
% % groups). You can study Indices, it just contains numbers between 1-10.
% % For each sample, we have a number which tells its group.
% Indices = crossvalind('Kfold', numel(Y), KFold);
% 
% % We have to run a loop where we study each of the groups (=KFold times)
% for ii=1:KFold
%     % The samples that do NOT belong to ii:th sample group are used when we
%     % build the PLS regression model in the ii:th run
%     ind_train = find(Indices ~= ii);
%     % The samples that do belong to ii:th sample group are not used when
%     % building the PLS regression model, but they are predicted, so they
%     % act as the validation samples
%     ind_val = find(Indices == ii);
% 
%     % First we build the ii:th model
%     [XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(X(ind_train,:),Y(ind_train),ncomp);
%     % Now we use this model to predict the samples that were left out from
%     % the training samples in this run. We store the predictions to
%     % yfit_PLS vector
%     yfit_PLS(ind_val) = [ones(numel(ind_val),1) X(ind_val,:)]*betaPLS;
% 
% end

[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,ncomp,'CV',10);
%yfit_PLS = [ones(numel(x_test),1) x_test]*betaPLS;
yfit_PLS_train = [ones(size(x_train,1),1) x_train]*betaPLS;
yfit_PLS_test = [ones(size(x_test,1),1) x_test]*betaPLS;
%yfit = [ones(size(x_test,1),1) x_test]*betaPLS;

% we can now plot a scatter plot and calculate the correlation between the
% reference values and predicted values. These results are now
% cross-validated, which means that we do have some idea about how the
% model actually performs when it is applied to data that has not been used
% when we build the model (well you know that we actually have used all the
% data when building the model, but still not all the data is used
% simultaneously)

figure
scatter(y_train,yfit_PLS_train)
title([' trainset r = ' num2str(corr(y_train,yfit_PLS_train))])

s1 = sum((y_train - yfit_PLS_train).^2);
s2 = sum((y_train - mean(yfit_PLS_train)).^2);
r2 = (1 - (s1/s2))*100
RMSEP = sqrt(mean((stats.Yresiduals).^2))
RMSECV = sqrt(msep(2,ncomp+1))

figure
scatter(y_test,yfit_PLS_test)
title([' testset r = ' num2str(corr(y_test,yfit_PLS_test))])

s11 = sum((y_test - yfit_PLS_test).^2);
s22 = sum((y_test - mean(yfit_PLS_test)).^2);
r22 = (1 - (s11/s22))*100

% we should also calculate the root-mean-squared error between y and
% yfit_PLS.

% Now that we have validated the model, we can build the "final model",
% where we use the whole X and Y data. This could be then applied to any
% new samples and we could say that the model performance is so and so
% (=typically meaning the pearson's correlation and rmsecv values)
% based on the earlier cross-validation.
[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(X,Y,ncomp);


%% Now we have calculated a cross-validated PLS regression model! 
% This sets the first baseline to our better models. We could
% repeat this same idea to all reference parameters that we have. One thing
% we have to take into account is that we do not have all reference
% parameters for all samples. So, for example, for equilibrium modulus, we
% need to find the samples for which we have the values:

X = Spectra_der2(find(~isnan(handles.equilibrium)),:);
Y = handles.equilibrium(find(~isnan(handles.equilibrium)));

% and then we could continue like shown above.