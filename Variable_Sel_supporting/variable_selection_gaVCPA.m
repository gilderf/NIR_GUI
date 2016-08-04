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

%% GA enhanced by VCP
%

A_max=10;
fold=5;
method='center';
Ratio_Better=0.1;
EDF_Run=50;
BMS_Run=1000;
A=5;


Result=VCPA(X,Y,A_max,fold,method,BMS_Run,EDF_Run,Ratio_Better);

selected_variables=Result.Vsel;

% then lets do a 10-fold cross-validation (we divide the data into 10
% groups) and create PLS model using only variables indicated by mcuve_vars
KFold = 10;
Indices = crossvalind('Kfold', numel(Y), KFold);

yfit_PLS_vcpa = zeros(size(Y));
for ii=1:KFold
    % these values are used for training a sub-model
    ind_train = find(Indices ~= ii);
    % these values are predicted
    ind_predict = find(Indices == ii);

    [XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(X(ind_train,selected_variables),Y(ind_train),A);
    % predicted values are stored
    yfit_PLS_vcpa(ind_predict) = [ones(numel(ind_predict),1) X(ind_predict,selected_variables)]*betaPLS;

end

% we can plot the result
figure
scatter(Y,yfit_PLS_vcpa)
title(num2str(corr(Y,yfit_PLS_vcpa)))