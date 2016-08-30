
figure
addpath supporting_files
clf;
load('Data_average_spectrums.mat');
load('Spectra_der2.mat');
load('UVE_complex.mat');
load('uve_selection.mat');

subplot(1,1,1) ; % main plot
set(gca,'DataAspectRatio',[1 1 1])
x1 = Wavelength;
y1 = Average_spectrums(90,:) ; %0.84 mm
y2 = Average_spectrums(70,:) ; %0.71 mm
y3 = Average_spectrums(184,:) ;%0.59 mm

plot(x1,y1,'-',x1,y2,'--',x1,y3,':','LineWidth',3);
set(gca, 'fontsize', 14, 'linewidth', 2)
set(gca, 'Ticklength', [0.009 0.0])
set(gca,'FontSize',12,'TickDir','out')

%grid on
lg = legend('0.84 mm','0.71 mm','0.59 mm','Location','northeast');
title(lg,'Thickness')
set(lg, 'EdgeColor', 'w')
title(lg,'Thickness')
xlim([700 1050])
ylim([0.3 1.3])
%ylim([0.2 1])
xlabel('Wavelength [nm]','fontWeight','bold')
% relative_offset = 1.5;
% xh = get(gca,'XLabel'); % Handle of the x label
% pause(0.2)
% set(xh, 'Units', 'Normalized')
% pause(0.2)
% pos = get(xh, 'Position');
% set(xh, 'Position',pos.*[1,relative_offset,1])
ylabel('Absorbance (A.U.)','fontWeight','bold')

%axes('position',[0.2 0.68 0.2 0.2]) ; % 1st inset
axes('position',[0.2 0.55 0.5 0.3]) ; % 1st inset
x1 = Wavelength(1:1998);
y1 = Spectra_der2(90,:) ;
y2 = Spectra_der2(70,:) ;
y3 = Spectra_der2(184,:) ;
plot(x1,y1,'-',x1,y2,'--',x1,y3,':','LineWidth',0.8);
%ylabel('Savits-Golay 2nd Der.','fontWeight','bold')
title('Savits-Golay 2nd Der.','fontWeight','bold')
grid on
%xlim([900 1000]);
xlim([920 980]);

axes('position',[0.13 0.1 0.775 0.11]) ; % 2nd inset
ind = wavelength_final +895;
s = zeros(size(Wavelength,2));
s(ind) = 1;
plot(Wavelength,0,'LineStyle','none','marker','none')
hold on
vline(ind,'k');
hold off
set(gca,'YTickLabel',{' '})
xlim([895 1540])
set(gca,'XTickLabel',{' '})
set(gca, 'Ticklength', [0 0])
%suptitle('Cartilage Thickness and UVE Selection Bands')
set(gcf, 'color', 'w')