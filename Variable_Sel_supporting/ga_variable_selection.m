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

combination = [X Y] ;

test1 = gaplsopt(combination,1) ;
test2 = gaplsopt(combination,2) ;

[b,c,d] = gaplssp(combination,182);

plotone(combination,b,c,d);
plotmore(combination,b);
plot(d);


%% TEST begins
A = 6;
[a,b] = sort(abs(d),'descend');

KFold = 10;
Indices = crossvalind('Kfold', numel(Y), KFold);

% here we have additional loop (jj) which adds variables one by one based
% on the reliability index values

% jj starts from A (number of PLS components), because that is the minimum
% number of variables we need when we use A PLS components.
for jj=A:size(X,2)
    for ii=1:KFold

        ind_predict = find(Indices == ii);
        ind_train = find(Indices ~= ii);

        [XL,yl,XS,YS,beta,PCTVAR] = plsregress(X(ind_train,b(1:jj)),Y(ind_train),A);
        yfit_GA_var(ind_predict,jj) = [ones(numel(ind_predict),1) X(ind_predict,b(1:jj))]*beta;

    end
end

%% Plot the correlation coefficient as a function of the number of variables
figure
plot(1:size(X,2),corr(yfit_GA_var,Y),'k','linewidth',2)
hold on
% this find the maximum correlation
[corr_coefficient, var_no] = max(corr(yfit_GA_var,Y));
% mark the best correlation
plot(var_no,corr_coefficient,'o')
hold off
xlabel('Number of variables','fontsize',18,'fontweight','bold')
ylabel('Correlation coefficient','fontsize',18,'fontweight','bold')
set(gca,'fontsize',16,'fontweight','bold')
xlim([0 size(X,2)])

%% we can now plot the best model because we know how many variables are needed

yfit = yfit_GA_var(:,var_no);

figure
scatter(Y,yfit)
title(num2str(corr(Y,yfit)))