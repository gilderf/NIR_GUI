
clear all;
close all;
clc;
%% Load and preprocess the data
% This files loads spectra and reference data, and does preprocessing to
% spectra.
% Needed files: data_average_spectrums.mat, Reference_data.xlsx
% NOTE: this is just copied from Jaakko's script. We may study also the
% preprocessing steps at some point.

% spectra are loaded, and wavelengths and spectra are stored in separate
% variables
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
%% BiPLS testing

% X = Spectra_der2;
xaxis = handles.wavelength(wave_start:wave_end);
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

bipls_var;

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
[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train(:,bipls_var),y_train,ncomp,'CV',10);
yfit_PLS_train = [ones(size(x_train,1),1) x_train(:,bipls_var)]*betaPLS;
yfit_PLS_test = [ones(size(x_test,1),1) x_test(:,bipls_var)]*betaPLS;


s1 = sum((y_train - yfit_PLS_train).^2);
s2 = sum((y_train - mean(yfit_PLS_train)).^2);
r2train(ncomp) = (1 - (s1/s2))*100;
rmsepC(ncomp) = sqrt(mean((stats.Yresiduals).^2));
rmsecvC(ncomp) = sqrt(msep(2,ncomp+1));

s11 = sum((y_test - yfit_PLS_test).^2);
s22 = sum((y_test - mean(yfit_PLS_test)).^2);
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
suptitle('biPLS')





