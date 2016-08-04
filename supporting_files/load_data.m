function [ X, Y_complex, wavelength] = load_data()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[fn,path] = uigetfile({'*.mat'}, 'Load spectroscopic data');
myData = importdata([path fn]);

% Setting Wavelength and Spectra as seperate matrixes, no error peaks
wavelength = myData.Wavelength(1:2000);
X = myData.Average_spectrums(:,1:2000);


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
% freq2 = myRef(:,9);
% freq3 = myRef(:,10);
% freq4 = myRef(:,11);
% freq5 = myRef(:,12);
% freq6 = myRef(:,13);
% freq7 = myRef(:,14);

% switch prop_selected{1}
%     case 'Cart. Thickness'
%         Y = thickness;
%     case 'Cart. Cal. Thickness'
%         Y = thickness_cal;
%     case 'Inst.Mod'
%         Y = instant;
%     case 'Dynamic.Mod'
%         Y = freq1;
%     case 'Equilibrium.Mod'
%         Y = equilibrium;
%     otherwise
%         warning('Invalid Property')
% end



end

