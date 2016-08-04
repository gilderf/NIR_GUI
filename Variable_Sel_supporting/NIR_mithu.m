%% Final NIR-analysis
% Getting the spectral and reference data
% [fn,path] = uigetfile({'*.mat'}, 'Load spectroscopic data');
myData = importdata('D:\NIR_GUI\Data_average_spectrums.mat');
handles.wavelength = myData.Wavelength(1:2000);
handles.spectra = myData.Average_spectrums(:,1:2000);
% [fn,path] = uigetfile({'*.xlsx'}, 'Load reference data');
myRef = xlsread('D:\NIR_GUI\Reference_data.xlsx');
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
ind_1 = 529; % Sample 10B - ICRS 0
ind_2 = 553;
ind_3 = 654; % Sample 10I - ICRS 1
ind_4 = 678;
ind_5 = 261; % Sample 5G - ICRS 2
ind_6 = 280;

% Components in analyses
% ncomp = 5;
% Creating visualizing images
%q = [180 356 532 708 884];
%Image = imread('D:\Biomechanics\Equine Official Indentation\Pohja_2.jpg');
%Image_base = zeros(size(Image,1));

% Font size
F_size = 12;
F_size_2 = 12;

%% Instant
ncomp = 6;
% Help variables
handles.spectra2 = Spectra_der2;
handles.spectra = handles.spectra2;
handles.ref = handles.instant./1e6;
% handles.ref = handles.thickness;
handles.ref2 = handles.ref;
minimi = min(handles.ref2);
maksimi = max(handles.ref2);
% Creating model and test set
handles.test = handles.spectra2([ind_1:ind_2 ind_3:ind_4 ind_5:ind_6],:);
handles.test_ref = handles.ref2([ind_1:ind_2 ind_3:ind_4 ind_5:ind_6],:);
handles.spectra([ind_1:ind_2 ind_3:ind_4 ind_5:ind_6],:) = [];
handles.ref([ind_1:ind_2 ind_3:ind_4 ind_5:ind_6],:) = [];
handles.model = handles.spectra;
handles.model_ref = handles.ref;
% Wavelength indexes
wave_start = 940;
wave_end = 1340;
%
% wave_start = 896;
% wave_end = 1540;

% Model setting
[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(handles.model(:,wave_start:wave_end),handles.model_ref,ncomp,'CV', 10);
respfitPLS = [ones(size(handles.model,1),1) handles.model(:,wave_start:wave_end)]*betaPLS;
ss_res = sum((handles.model_ref - respfitPLS).^2);
ss_tot = sum((handles.model_ref - mean(handles.model_ref)).^2);
r2 = (1 - (ss_res/ss_tot))*100;
RMSEP = sqrt(mean((stats.Yresiduals).^2));
RMSECV = sqrt(msep(2,ncomp+1));
% Test setting
testpred = [ones(size(handles.test,1),1) handles.test(:,wave_start:wave_end)]*betaPLS;
ss_res_2 = sum((handles.test_ref - testpred).^2);
ss_tot_2 = sum((handles.test_ref - mean(handles.test_ref)).^2);
r22 = (1 - (ss_res_2/ss_tot_2))*100;
RMSEP_2 = sqrt(mean((handles.test_ref-testpred).^2));


% Plotting model and test
subplot(2,1,1); plot(handles.model_ref,respfitPLS,'o','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
set(gca,'FontSize',F_size,'fontWeight','bold')
xlim([0 25]); ylim([0 20])
xlabel('Measured inst. modulus [MPa]','FontSize',F_size,'fontWeight','bold')
ylabel('Predicted inst. modulus [MPa]','FontSize',F_size_2,'fontWeight','bold')
h = legend(sprintf('R^{2} = %5.2f, N = %3.0f, p < 0.001',round(100*r2)./100,size(handles.model,1))); 
h2 = get(h,'children'); 
delete(h2(1)) ;
legend(gca,'boxoff')
ylab = get(gca,'YLabel');
handles2(5).ylab = get(ylab,'Position');
set(ylab,'Position',get(ylab,'Position')-[0.28 0 0])
p = polyfit(handles.model_ref,respfitPLS,1);
y = p(1).*[0 20]+p(2);
hold on
plot([0 25],y,'--k','LineWidth',2)

subplot(2,1,2); plot(handles.test_ref,testpred,'o','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
set(gca,'FontSize',F_size,'fontWeight','bold')
xlim([0 25]); ylim([0 20])
xlabel('Measured inst. modulus [MPa]','FontSize',F_size,'fontWeight','bold')
ylabel('Predicted inst. modulus [MPa]','FontSize',F_size_2,'fontWeight','bold')
ylab = get(gca,'YLabel');
handles2(6).ylab = get(ylab,'Position');
set(ylab,'Position',get(ylab,'Position')-[0.28 0 0])
p = polyfit(handles.test_ref,testpred,1);
y = p(1).*[0 20]+p(2);
hold on
plot([0 25],y,'--k','LineWidth',2)
h = legend(sprintf('R^{2} = %5.2f, N = %3.0f, p < 0.001',round(100*r22)./100,size(handles.test,1))); 
h2 = get(h,'children'); 
delete(h2(1));
legend(gca,'boxoff')