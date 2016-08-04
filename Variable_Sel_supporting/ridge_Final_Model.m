%Lasso Model building
ccc;
%% Load and preprocess the data
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
% The best test set was found to be Sample 10B, Sample 10I and Sample 5G
% ind_1 = 529; % Sample 10B - ICRS 0
% ind_2 = 553;
% ind_3 = 654; % Sample 10I - ICRS 1
% ind_4 = 678;
% ind_5 = 261; % Sample 5G - ICRS 2
% ind_6 = 280;

% Wavelength indexes
%Only for Instant Modulus
% wave_start = 940;
% wave_end = 1340;

%For rest of parameters
wave_start = 896;
wave_end = 1540;

%% Seperating the data into training and testing data.
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


% x_train = x_train(1:10,:);
% y_train = y_train(1:10,:);
%
% x_test = x_test(1:10,:);
% y_test = y_test(1:10,:);

% Perform Ridge regression

% Find optimum lambda
interim = [0 : .1 : 1];
till = size(interim,2);

r2train = zeros(till,1);
RMSCV = zeros(till,1);
r2test = zeros(till,1);
RMSEP = zeros(till,1);

C = zeros(till,1);

for i = 1:till
    
    [B] = ridge(y_train,x_train,interim(i),0);
    
    y_predicted_train = B(1) + x_train*B(2:end);
    y_predicted_test = B(1) + x_test*B(2:end);
    
    %corr(y_predicted,y_train)
    %C(i)= corr(y_predicted_test,y_test)
    
    s1 = sum((y_train - y_predicted_train).^2);
    s2 = sum((y_train - mean(y_predicted_train)).^2);
    r2train(i) = (1 - (s1/s2))*100;
    RMSCV(i) = sqrt(mean((y_train-y_predicted_train).^2));
    
    s11 = sum((y_test - y_predicted_test).^2);
    s22 = sum((y_test - mean(y_predicted_test)).^2);
    r2test(i) = (1 - (s11/s22))*100;
    RMSEP(i) = sqrt(mean((y_test-y_predicted_test).^2));
    
    
    
    
end

% Print the results Grid
ncompmax = till;
figure;
title('Baseline per manuscript');
subplot(2,2,1);
plot(1:ncompmax,r2train)
xlabel('Lambda Value');
ylabel('R^2 Train');
hold on
plot(find(r2train == min(r2train)),min(r2train),'o')
plot(find(r2test == max(r2test)),r2train,'X')
hold off

subplot(2,2,3);
plot(1:ncompmax,r2test)
xlabel('Lambda Value');
ylabel('R^2 Test');
hold on
plot(find(r2test == max(r2test)),max(r2test),'o')
hold off

subplot(2,2,4);
plot(1:ncompmax,RMSEP)
xlabel('Lambda Value');
ylabel('RMSEP');
hold on
plot(find(RMSEP == min(RMSEP)),min(RMSEP),'o')
plot(find(r2test == max(r2test)),RMSEP,'X')
hold off

subplot(2,2,2);
plot(1:ncompmax,RMSCV)
xlabel('Lambda Value');
ylabel('RMSECV');
hold on
plot(find(RMSCV == min(RMSCV)),min(RMSCV),'o')
plot(find(r2test == max(r2test)),RMSCV,'X')
hold off
suptitle('Ridge Regression')

