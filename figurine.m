function [  ] = figurine( func_prop_sel )
% Plots figure for MC-UVE-PLS variable selection method
%   Detailed explanation goes here
%ccc;

load('Data_average_spectrums.mat');
load('Spectra_der22.mat');
load('UVE_complex.mat');
load('uve_selected_var.mat');
%load('Range.mat');
%'wave_start','wave_end'
%close all;
%clc;
s = abs(UVE.RI);
s1 = s - 1.24 ;
ind = find(s1 < 0);
s1 = s1+1.2;
s1 (ind) = 0;
figure
subplot(5,1,5)
s = abs(UVE.RI);
s1 = s - 1.24 ;
ind = find(s1 < 0);
s1 = s1+1.2;
s1 (ind) = 0;
ind2 = find(s1 ~= 0);
s1(ind2) = 1;
%axis off
grid on
bar((896 : 1540),s1);
%bar(s1)
hold on
%ylabel('Variable Importance','fontWeight','bold')
xlabel('MC-UVE-PLS selection band','fontWeight','bold')
set(gca,'YTickLabel',{' '})
subplot(5,1,[3 4])
x1 = ( 800 : 1600 );
y1 = Spectra_der22(90,800 : 1600) ;
y2 = Spectra_der22(70,800 : 1600) ;
y3 = Spectra_der22(184,800 : 1600) ;
%plot((800 : 1600),Spectra_der22(106,800:1600),'DisplayName','Spectra_der22')
plot(x1,y1,x1,y2,x1,y3,'LineWidth',0.8);
ylabel('Double Savits-Golay 2nd Der.','fontWeight','bold')
grid on
% hold on
% bar((896 : 1540),not(s1));
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','r','facealpha',0.15,'EdgeAlpha',0.15)
% hold off
subplot(5,1,[1 2])
%plot((896 : 1521),Average_spectrums(106,896:1521),'DisplayName','Average_spectrums');
%plot((896 : 1521),Average_spectrums(106,896:1521),'color',[1 1 1]);
x1 = ( 896 : 1521 );
y1 = Average_spectrums(90,896:1521) ; %0.84 mm
y2 = Average_spectrums(70,896:1521) ; %0.71 mm
y3 = Average_spectrums(184,896:1521) ;%0.59 mm
%plot((896 : 1521),Average_spectrums(106,896:1521),'b','LineWidth',3);
%plot(x1,y1,'b','LineWidth',3,x1,y2,'LineWidth',3);
plot(x1,y1,x1,y2,x1,y3,'LineWidth',3);
%plot((896 : 1521),Average_spectrums(106,896:1521),'w');
hold on
bar((896 : 1540),not(s1));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','r','facealpha',0.15,'EdgeAlpha',0.15)
hold off
grid on
xlabel('Absorbance NIR Values','fontWeight','bold')
%hold on
%plot((find(s1 ~= 0)+895),Spectra_der2(106,(find(s1 ~= 0)+895)),'X')
ylabel('Absorbance NIR Values','fontWeight','bold')
xlabel('Wavelength [nm]','fontWeight','bold')
title_to_display = strcat('MC-UVE-PLS : ',func_prop_sel);
%suptitle(title_to_display)
suptitle(title_to_display)

end

