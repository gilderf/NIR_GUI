%plotmore plots the variable selections from different programs
%run on the same data set
%It requires the data set and a matrix created by adding the b vectors
%from the programs: tot=[b1;b2;...;bn] 
function plotmore(dataset,tot)
[r,co]=size(tot);
dataset=dataset(:,1:co);
mi=min(min(dataset));
ma=max(max(dataset));
figure
plot(dataset','r')
hold on
set(gca,'XLim',[.5 co+.5]);
ra=ma-mi;
set(gca,'YLim',[mi-ra/20 ma+ra/20])
for i=1:r
  lev=mi-(ra/20)*i;  disp([' '])
  disp([' '])
  disp(['Program ' int2str(i)])
  % Code added by Mithilesh Prakash on 2,8,2016 for quickening the GA process
  %re=input('How many variables? ');
  re= 100;
  for ii=1:re
    plot([tot(i,ii)-.5 tot(i,ii)+.5],[lev lev],'b');
  end
end
set(gca,'YLim',[lev-ra/20 ma+ra/20])
hold off
figure(gcf)
zoomrb


