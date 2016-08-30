

% UVE.RI shows us the reliability of the variables.
figure
plot(abs(UVE.RI),'linewidth',2);
xlabel('variable index');
ylabel('reliability index');
set(gcf,'color','w');
hold on
[p,locs] = findpeaks(abs(UVE.RI));
plot(locs(find(p > 1.24)),p(find(p > 1.24)),'ro');
hold on
plot(x_train(find(p > 1.24)),'DisplayName','x_train')

%****
close all ;
clc;
s = abs(UVE.RI);
s1 = s - 1.24 ;
ind = find(s1 < 0);
s1 = s1+1.2;
s1 (ind) = 0;
ind2 = find(s1 ~= 0);
s1(ind2) = 1;
figure
bar((896 : 1540),s1);
%bar(s1)
hold on
ylabel('Variable Importance','fontWeight','bold')
xlabel('Wavelength','fontWeight','bold')


%***
%for i = 1 : size(Spectra_der2,1)
close all;
%figure
subplot(3,1,1)
plot(Spectra_der2(1,:),'DisplayName','Spectra_der2')
subplot(3,1,2)
plot((800 : 1600),Spectra_der2(106,800:1600),'DisplayName','Spectra_der2')
subplot(3,1,3)
s = abs(UVE.RI);
s1 = s - 1.24 ;
ind = find(s1 < 0);
s1 = s1+1.2;
s1 (ind) = 0;
%figure
bar((896 : 1540),s1);
hold on
%k = waitforbuttonpress;
%i
%end

close all;
clc;
s = abs(UVE.RI);
s1 = s - 1.24 ;
ind = find(s1 < 0);
s1 = s1+1.2;
s1 (ind) = 0;
%figure
subplot(3,1,1)
plot((800 : 1600),Spectra_der2(106,800:1600),'DisplayName','Spectra_der2')
grid on
hold on
plot((find(s1 == 0)+895),Spectra_der2(106,(find(s1 == 0)+895)),'X')
ylabel('MeasuredM inst. modulus [MPa]','fontWeight','bold')
xlabel('Wavelength [nm]','fontWeight','bold')
subplot(3,1,2)
%plot((800 : 1600),Spectra_der2(106,800:1600),'DisplayName','Spectra_der2')
grid on
%hold on
plot((find(s1 ~= 0)+895),Spectra_der2(106,(find(s1 ~= 0)+895)),'X')
ylabel('MeasuredM inst. modulus [MPa]','fontWeight','bold')
xlabel('Wavelength [nm]','fontWeight','bold')
subplot(3,1,3)
%figure
bar((896 : 1540),s1);
grid on
ylabel('Variable Importance Level')
% ax = gca;
% ax.XTick = 800 : 1600 ;
% ax.XTickLabels = 800 : 1600;
%
% legend({'Average High','Average Low'},'Location','northwest')

%% New
close all;
clc;
s = abs(UVE.RI);
s1 = s - 1.24 ;
ind = find(s1 < 0);
s1 = s1+1.2;
s1 (ind) = 0;
%figure
subplot(4,1,4)
s = abs(UVE.RI);
s1 = s - 1.24 ;
ind = find(s1 < 0);
s1 = s1+1.2;
s1 (ind) = 0;
ind2 = find(s1 ~= 0);
s1(ind2) = 1;
grid on
bar((896 : 1540),s1);
%bar(s1)
hold on
ylabel('Variable Importance','fontWeight','bold')
xlabel('Wavelength','fontWeight','bold')
subplot(4,1,[1 3])
plot((800 : 1600),Spectra_der2(106,800:1600),'DisplayName','Spectra_der2')
grid on
%hold on
%plot((find(s1 ~= 0)+895),Spectra_der2(106,(find(s1 ~= 0)+895)),'X')
ylabel('MeasuredM inst. modulus [MPa]','fontWeight','bold')
xlabel('Wavelength [nm]','fontWeight','bold')
% subplot(3,1,3)
% s = abs(UVE.RI);
% s1 = s - 1.24 ;
% ind = find(s1 < 0);
% s1 = s1+1.2;
% s1 (ind) = 0;
% %figure
% bar((896 : 1540),s1);
% grid on
% ylabel('Variable Importance Level')

%% triple section
close all;
clc;
s = abs(UVE.RI);
s1 = s - 1.24 ;
ind = find(s1 < 0);
s1 = s1+1.2;
s1 (ind) = 0;
%figure
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
plot((800 : 1600),Spectra_der2(106,800:1600),'DisplayName','Spectra_der2')
ylabel('Savits-Golay 2nd Der.','fontWeight','bold')
grid on
subplot(5,1,[1 2])
plot((896 : 1521),Average_spectrums(106,896:1521),'DisplayName','Average_spectrums');
grid on
xlabel('Absorbance NIR Values','fontWeight','bold')
%hold on
%plot((find(s1 ~= 0)+895),Spectra_der2(106,(find(s1 ~= 0)+895)),'X')
ylabel('Absorbance NIR Values','fontWeight','bold')
xlabel('Wavelength [nm]','fontWeight','bold')
