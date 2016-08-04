function ipcascoplotall(Model,maxPC,samplenames,classvariable,no_of_plots_pr_figure)

%  ipcascoplotall makes all 2D score plots up to maxPC for all intervals
%     NOTE: Press [enter] to continue to next combination
%
%  Input:
%  Model (the output from ipca.m)
%  maxPC: up to this number all combination of score plots are produced
%  samplenames: string, see makeSampleNames for an example (use [] if not available)
%  classvariable: vector of numbers (1 to n) assigning each class (1 to n) 
%                 see makeClasses for an example (use [] if not available)
%                 up to n=7 classes can be coloured
%  no_of_plot_pr_figure must be 2, 4, 6 or 8 (optional, default is 6)
%
%  Copyright, Chemometrics Group - KVL, Copenhagen, Denmark
%  Lars N�rgaard, July 2004
%
%  ipcascoplotall(Model,maxPC,samplenames,classvariable,no_of_plots_pr_figure);

if nargin==0
   disp(' ')
   disp(' ipcascoplotall(Model,maxPC,samplenames,classvariable,no_of_plots_pr_figure);')
   disp(' ')
   disp(' Example:')
   disp(' ipcascoplotall(Model,4,[],classes);')
   disp(' ')
   return
end

if nargin<5
    no_of_plots_pr_figure=6;
end

set(0,'Units','pixels');
Scrsiz=get(0,'ScreenSize');
ScrLength=Scrsiz(3);
ScrHight=Scrsiz(4);
bdwidth=10;
% [left(->) bottom(up) width hight]
figpos=[0.44*ScrLength 0.05*ScrHight 0.16*ScrLength 0.05*ScrHight];

for i=1:maxPC
  for j=(i+1):maxPC
    iscoplot(Model,i,j,samplenames,classvariable,no_of_plots_pr_figure)
    h=msgbox('Press [Enter] when you want to continue with next combination','Continue','replace');
    set(h,'Position',figpos)
    %keyboard
    %get(h)
    pause
  end
end
close all
