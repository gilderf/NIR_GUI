function [ xtr, ytr, xtt, ytt, range ] = preprocess_nir_data( X, Y, wavelength,continue_flag, func_prprty )
%UNTITLED2 Summary of this function goes here
% This does preprocessing to spectra.
% First, Savitzky-Golay smoothing is applied. If you are not yet familiar
% with Savitzky-Golay, you should read a little bit about it. It is quite
% simple, and a very common method in spectroscopy.

% Second, the second derivative spectra area calculated

% Third, a moving-average filter is applied to reduce noise.

% Best smoothing    1) window 40, degree 3
%                   2) Second derivative
%                   3) window 10, degree 5



if (continue_flag == 0)
    % h = waitbar(0,'Please wait...','Name','Preprocessing NIR Data...');
    h = waitbar(0);
    for kk = 1:size(X,1)
        Spectra_smooth = smooth(wavelength,X(kk,:),40,'sgolay',3);
        Spectra_smooth = Spectra_smooth';
        for ll = 1:size(Spectra_smooth,2)-2
            Spectra_der(:,ll) = (Spectra_smooth(:,ll+2)-2*Spectra_smooth(:,ll+1)+Spectra_smooth(:,ll))./(((wavelength(ll+2)-wavelength(ll))/2).^2);
        end
        wind = 10; % Normal 10
        Spectra_der2(kk,:) = filter(ones(1,wind)/wind,1,Spectra_der);
        waitbar(kk / size(X,1));
    end
     save('Spectra_der2.mat','Spectra_der2');
   % Extra processing for double Savits-Golay, for figure plotting in
   % paper1, remove while analysing or add option in GUI
        for kk = 1:size(Spectra_der2,1)
            Spectra_smooth = smooth(wavelength(1,1:1998),Spectra_der2(kk,:),40,'sgolay',3);
            Spectra_smooth = Spectra_smooth';
            for ll = 1:size(Spectra_smooth,2)-2
                Spectra_der(:,ll) = (Spectra_smooth(:,ll+2)-2*Spectra_smooth(:,ll+1)+Spectra_smooth(:,ll))./(((wavelength(ll+2)-wavelength(ll))/2).^2);
            end
            wind = 10; % Normal 10
            Spectra_der22(kk,:) = filter(ones(1,wind)/wind,1,Spectra_der);
            %waitbar(kk / size(Spectra_der2,1));
        end
    %end
    close(h)
    clear Spectra_smooth kk ll  der Spectra_der
   
    %only for Diagram not used in spectral analysis
    save('Spectra_der22.mat','Spectra_der22');
end

if (continue_flag == 1 && isequal(Y,0))
    [fn,path] = uigetfile({'*.xlsx'}, 'Load reference data');
    %myRef = xlsread('D:\NIR_GUI\Reference_data.xlsx');
    myRef = xlsread([path fn]);
    Y_complex.thickness = myRef(:,1);
    Y_complex.thickness_cal = myRef(:,2);
    Y_complex.ICRS = myRef(:,3);
    Y_complex.instant = myRef(:,5);
    Y_complex.instant = Y_complex.instant./1e6;
    Y_complex.equilibrium = myRef(:,7);
    Y_complex.equilibrium = Y_complex.equilibrium./1e6;
    Y_complex.freq1 = myRef(:,8);
    Y_complex.freq1 = Y_complex.freq1./1e6;
    [ Y ] = property_selector( Y_complex,func_prprty );
end
% The best test set was found to be Sample 10B, Sample 10I and Sample 5G
% ind_1 = 529; % Sample 10B - ICRS 0
% ind_2 = 553;
% ind_3 = 654; % Sample 10I - ICRS 1
% ind_4 = 678;
% ind_5 = 261; % Sample 5G - ICRS 2
% ind_6 = 280;

% Wavelength indexes
wave_start = 896;
wave_end = 1540;

if ( isequal(func_prprty{1},'Inst.Mod') == 1)
    %Only for Instant Modulus
    wave_start = 940;
    wave_end = 1340;
end
if ( isequal(func_prprty{1},'Equilibrium.Mod') == 1)
    %for equilibrium modulus
    wave_start = 900;
    wave_end = 1524;
end

%Uncomment for full spectrum
%     wave_start = 1;
%     wave_end = 1998;
range = wave_start : wave_end ;

% if exist('Range','file') == 2
% delete('Range.mat');
% end
% save('Range.mat','wave_start','wave_end');

%% Seperating the data into training and testing data.
% In this case, lets use second derivative spectra (we could use the
% original spectra without preprocessing as well)
% Spectra_der2 = load('Spectra_der2.mat');
% X = Spectra_der2.Spectra_der2(:,wave_start:wave_end);
Spectra_der2 = load('Spectra_der2.mat');
X = Spectra_der2.Spectra_der2(:,wave_start:wave_end);
x_train =[X(1:260,:) ; X(281:528,:) ; X(554:653,:) ; X(679:869,:)];
x_test = [X(261:280,:) ;X(529:553,:) ;X(654:678,:)];

% And the predicted feature is instant modulus
% Y = handles.instant;
% Y = Y./1e6;
%Y = handles.equilibrium;
y_train =[Y(1:260,:) ; Y(281:528,:) ; Y(554:653,:) ; Y(679:869,:)];
y_test = [Y(261:280,:) ;Y(529:553,:) ;Y(654:678,:)];

xtr = x_train(find(~isnan(y_train)),:);
ytr = y_train(find(~isnan(y_train)));

xtt = x_test(find(~isnan(y_test)),:);
ytt = y_test(find(~isnan(y_test)));

%Adding ZScore processing

% xtr = zscore(xtr);
% ytr = zscore(ytr);
% xtt = zscore(xtt);
% ytt = zscore(ytt);


end

