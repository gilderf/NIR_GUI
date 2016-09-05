figure
addpath supporting_files
clf;
load('D:\NIR Gui Project\Data_average_spectrums.mat');
load('D:\NIR Gui Project\Spectra_der2.mat');
load('D:\NIR Gui Project\UVE_complex.mat');
load('D:\NIR Gui Project\uve_selection.mat');

subplot(1,1,1) ; % main plot
set(gca,'DataAspectRatio',[1 1 1])
x1 = Wavelength;
y1 = Average_spectrums(97,:) ; %2.20 MPa 
y2 = Average_spectrums(68,:) ; %1.58 MPa
y3 = Average_spectrums(191,:) ;%1.02 MPa 

x1 = x1(1:10:end);
y1 = y1(1:10:end);
y2 = y2(1:10:end);
y3 = y3(1:10:end);
plot(x1,y1,'r','LineStyle','-','LineWidth',3,'MarkerSize',2,'marker','none');
hold on
plot(x1,y2,'b','LineStyle','--','LineWidth',3,'MarkerSize',2,'marker','none');
hold on
plot(x1,y3,'k','LineStyle','-','LineWidth',1,'MarkerSize',4,'marker','*');
set(gca, 'fontsize', 14, 'linewidth', 2)
set(gca, 'Ticklength', [0.009 0.0])
set(gca,'FontSize',12,'TickDir','out','fontWeight','bold')

%grid on
lg = legend('2.20 MPa','1.58 MPa','1.02 MPa','Location','northeast','fontsize', 14);
title(lg,'Thickness')
set(lg, 'EdgeColor', 'w')
title(lg,'Thickness')
xlim([700 1050])
ylim([0.25 1.25])
xlabel('Wavelength (nm)','fontWeight','bold','fontsize', 16)
ylabel('Absorbance (A.U.)','fontWeight','bold','fontsize', 16)
axes('position',[0.2 0.55 0.46 0.3]) ; % 1st inset
x1 = Wavelength(1:1998);
y1 = Spectra_der2(97,:) ;
y2 = Spectra_der2(68,:) ;
y3 = Spectra_der2(191,:) ;
plot(x1,y1,'r','LineStyle','-','LineWidth',1,'MarkerSize',2,'marker','none');
hold on
plot(x1,y2,'b','LineStyle','--','LineWidth',2,'MarkerSize',2,'marker','none');
hold on
plot(x1,y3,'k','LineStyle','-','LineWidth',1,'MarkerSize',2,'marker','*');
title('Savits-Golay 2nd Der.','fontWeight','bold','fontsize', 14)
set(gca,'FontSize',12,'TickDir','out','fontWeight','bold')
xlim([925 975]);

axes('position',[0.13 0.1 0.775 0.11]) ; % 2nd inset
ind = wavelength_final +899;
s = zeros(size(Wavelength,2));
s(ind) = 1;
plot(Wavelength,0,'LineStyle','none','marker','none')
hold on
vline(ind,'k');
hold off
set(gca,'YTickLabel',{' '})
xlim([899 1523])
set(gca,'XTickLabel',{' '})
set(gca, 'Ticklength', [0 0])
set(gcf, 'color', 'w')
%suptitle('E')
%%
%set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPositionMode','auto')
% print('test', '-dtiff', '-r300');