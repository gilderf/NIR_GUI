function plspvsm(Model,no_of_int_lv,interval,y_variable)

%  plspvsm plots predicted versus measured for a combination of several intervals
%
%  Input:
%  Model is the output from ipls.m, plsmodel.m or plspredict.m
%  no_of_int_lv is the number of PLS components to use for the interval model
%  interval: should be given if ipls model is input, otherwise state [] or
%    omit. Use 0 for global model.
%  y_variable is the number of the y-variable that the plot is made for
%    in the case of only one y-variable simply omit or type 1
%
%  Copyright, Chemometrics Group - KVL, Copenhagen, Denmark
%  Lars N�rgaard, July 2004
%
%  plspvsm(Model,no_of_int_lv);

if nargin==0
   disp(' ')
   disp(' plspvsm(Model,no_of_int_lv,interval,y_variable);')
   disp(' ')
   disp(' Example:')
   disp(' plspvsm(Model,5,10,1);')
   disp(' ')
   disp(' plspvsm(Model,5);')
   disp(' ')
   return
end

if ~ismember(Model.type,{'PLS','iPLS','PLSprediction'})
    disp(' ')
    disp('This function only works with output from ipls.m, plsmodel.m or plspredict.m')
    disp(' ')
    return
end

if strcmp(Model.type,'iPLS') & nargin<=2
    disp(' ')
    disp('Plotting results from iPLS model: Remember to give interval number as the third parameter')
    disp(' ')
    return
end

if nargin>=3
  if ismember(Model.type,{'PLS','PLSprediction'}) & ~isempty(interval)
     disp(' ')
     disp('Plotting results from PLS/PLSprediction model: It is not necessary to specify interval')
     disp('Use [] or omit if last parameter')
     disp(' ')
  end
end

if nargin<=3
   y_variable=1;
end

if nargin >=3 
  if interval==0 & strcmp(Model.type,'iPLS')
     interval=Model.intervals+1;
  elseif interval==0 & strcmp(Model.type,'PLSprediction')
     interval=Model.CalModel.intervals+1;
  end
end

if nargin==2 & strcmp(Model.type,'PLSprediction')
   interval=Model.CalModel.intervals+1;
end

set(0,'Units','pixels');
Scrsiz=get(0,'ScreenSize');
ScrLength=Scrsiz(3);
ScrHight=Scrsiz(4);
bdwidth=10;
% [left(->) bottom(up) width hight]
pos1=[bdwidth (0.4*ScrHight+bdwidth) (ScrLength/2-2*bdwidth) ScrHight/1.7-(70+bdwidth)];
pos2=[pos1(1)+ScrLength/2 pos1(2) pos1(3) pos1(4)];

% Position of interval(s)
figure('Position',pos1);

switch Model.type
  case 'iPLS'
	if isempty(Model.xaxislabels)
       plot(Model.rawX','k')
       xlabel('Variables')
       stvar=Model.allint(interval,2);
       endvar=Model.allint(interval,3);
       if interval<size(Model.allint,1)
          titletext=sprintf('Interval number %g, variables %g-%g',interval,stvar,endvar);
       else
          titletext=sprintf('Global model, variables %g-%g',stvar,endvar);
       end
	else
       plot(Model.xaxislabels,Model.rawX,'k')
       xlabel('Wavelength')
       stwav=Model.xaxislabels(Model.allint(interval,2));
       endwav=Model.xaxislabels(Model.allint(interval,3));
       if interval<size(Model.allint,1)
          titletext=sprintf('Interval number %g, wavelengths %g-%g',interval,stwav,endwav);
       else
          titletext=sprintf('Global model, wavelengths %g-%g',stwav,endwav);
       end
	end
	ytext=sprintf('Response, raw data [%s is used in the calculations]', Model.prepro_method);
    ylabel(ytext);
	
	title(titletext)
	axis tight
	actualaxis=axis;
	hold on
	a=Model.allint(interval,2);
	b=Model.allint(interval,3);
	if isempty(Model.xaxislabels)
       h1temp=area([a b],[actualaxis(3) actualaxis(3)]); % Negative areas
       set(h1temp,'FaceColor',[0.75 0.75 0.75]);
       h2temp=area([a b],[actualaxis(4) actualaxis(4)]);
       set(h2temp,'FaceColor',[0.75 0.75 0.75]);
       plot(Model.rawX','k') % To overlay spectra on area plot
	else
       h1temp=area([Model.xaxislabels(a) Model.xaxislabels(b)],[actualaxis(3) actualaxis(3)]); % Negative areas
       set(h1temp,'FaceColor',[0.75 0.75 0.75]);
       h2temp=area([Model.xaxislabels(a) Model.xaxislabels(b)],[actualaxis(4) actualaxis(4)]);
       set(h2temp,'FaceColor',[0.75 0.75 0.75]);
       plot(Model.xaxislabels,Model.rawX,'k') % To overlay spectra on area plot
	end
	hold off
      
  case {'PLS'}
	if isempty(Model.xaxislabels)
       plot(Model.rawX','k')
       xlabel('Variables')
       if isfield(Model,'windowsize') % If oneModel is based on mwModel
           stvar=Model.selected_vars(1);
           endvar=Model.selected_vars(end);
           titletext=sprintf('Selected variables [%g to %g]',Model.selected_vars(1),Model.selected_vars(end));
       else
           stvar=Model.allint(Model.selected_intervals,2);
           endvar=Model.allint(Model.selected_intervals,3);
           titletext=sprintf('Selected intervals [%s]',num2str(Model.selected_intervals));           
       end
    else
       plot(Model.xaxislabels,Model.rawX,'k')
       xlabel('Wavelength')
       if isfield(Model,'windowsize') % If oneModel is based on mwModel
           stwav=Model.xaxislabels(Model.selected_vars(1));
           endwav=Model.xaxislabels(Model.selected_vars(end));
           titletext=sprintf('Selected variables [%g to %g]',Model.selected_vars(1),Model.selected_vars(end));
       else
           stwav=Model.xaxislabels(Model.allint(Model.selected_intervals,2));
           endwav=Model.xaxislabels(Model.allint(Model.selected_intervals,3));
           titletext=sprintf('Selected intervals [%s]',num2str(Model.selected_intervals));
       end       
    end
	ytext=sprintf('Response, raw data [%s is used in the calculations]', Model.prepro_method);
    ylabel(ytext);
	title(titletext)
	axis tight
	actualaxis=axis;
	hold on
	if isempty(Model.xaxislabels)
       if isfield(Model,'windowsize') % If oneModel is based on mwModel
           a=Model.selected_vars(1);
           b=Model.selected_vars(end);
       else
           a=Model.allint(Model.selected_intervals,2);
           b=Model.allint(Model.selected_intervals,3);
       end
       for i=1:max(size(a))
          h1temp=area([a(i) b(i)],[actualaxis(3) actualaxis(3)]); % Negative areas
          set(h1temp,'FaceColor',[0.75 0.75 0.75]);
          h2temp=area([a(i) b(i)],[actualaxis(4) actualaxis(4)]);
          set(h2temp,'FaceColor',[0.75 0.75 0.75]);
       end
       plot(Model.rawX','k') % To overlay spectra on area plot
	else
       if isfield(Model,'windowsize') % If oneModel is based on mwModel
           a=Model.xaxislabels(Model.selected_vars(1));
           b=Model.xaxislabels(Model.selected_vars(end));
       else
           a=Model.xaxislabels(Model.allint(Model.selected_intervals,2));
           b=Model.xaxislabels(Model.allint(Model.selected_intervals,3));
       end
       for i=1:max(size(a))
          h1temp=area([a(i) b(i)],[actualaxis(3) actualaxis(3)]); % Negative areas
          set(h1temp,'FaceColor',[0.75 0.75 0.75]);
          h2temp=area([a(i) b(i)],[actualaxis(4) actualaxis(4)]);
          set(h2temp,'FaceColor',[0.75 0.75 0.75]);
       end
		plot(Model.xaxislabels,Model.rawX,'k') % To overlay spectra on area plot
	end
	hold off
    
  case {'PLSprediction'}
	if isempty(Model.CalModel.xaxislabels)
       plot(Model.CalModel.rawX','k')
       xlabel('Variables')
       if isfield(Model.CalModel,'windowsize') % If predModel/oneModel is based on mwModel
           stvar=Model.CalModel.selected_vars(1);
           endvar=Model.CalModel.selected_vars(end);
           titletext=sprintf('Selected variables [%g to %g]',Model.CalModel.selected_vars(1),Model.CalModel.selected_vars(end));           
       else
           stvar=Model.CalModel.allint(Model.CalModel.selected_intervals,2);
           endvar=Model.CalModel.allint(Model.CalModel.selected_intervals,3);
           titletext=sprintf('Selected intervals [%s]',num2str(Model.CalModel.selected_intervals));
       end
	else
       plot(Model.CalModel.xaxislabels,Model.CalModel.rawX,'k')
       xlabel('Wavelength')
       if isfield(Model.CalModel,'windowsize') % If oneModel is based on mwModel
           stwav=Model.CalModel.xaxislabels(Model.CalModel.selected_vars(1));
           endwav=Model.CalModel.xaxislabels(Model.CalModel.selected_vars(end));
           titletext=sprintf('Selected variables [%g to %g]',Model.CalModel.selected_vars(1),Model.CalModel.selected_vars(end));
       else
           stwav=Model.CalModel.xaxislabels(Model.CalModel.allint(Model.CalModel.selected_intervals,2));
           endwav=Model.CalModel.xaxislabels(Model.CalModel.allint(Model.CalModel.selected_intervals,3));
           titletext=sprintf('Selected intervals [%s]',num2str(Model.CalModel.selected_intervals));
       end
	end
	ytext=sprintf('Response, raw data [%s is used in the calculations]', Model.CalModel.prepro_method);
    ylabel(ytext);
	title(titletext)
	axis tight
	actualaxis=axis;
	hold on
	if isempty(Model.CalModel.xaxislabels)
       if isfield(Model.CalModel,'windowsize') % If oneModel is based on mwModel
           a=Model.CalModel.selected_vars(1);
           b=Model.CalModel.selected_vars(end);
       else
           a=Model.CalModel.allint(Model.CalModel.selected_intervals,2);
           b=Model.CalModel.allint(Model.CalModel.selected_intervals,3);
       end
       for i=1:max(size(a))
          h1temp=area([a(i) b(i)],[actualaxis(3) actualaxis(3)]); % Negative areas
          set(h1temp,'FaceColor',[0.75 0.75 0.75]);
          h2temp=area([a(i) b(i)],[actualaxis(4) actualaxis(4)]);
          set(h2temp,'FaceColor',[0.75 0.75 0.75]);
       end
       plot(Model.CalModel.rawX','k') % To overlay spectra on area plot
	else
       if isfield(Model.CalModel,'windowsize') % If oneModel is based on mwModel
           a=Model.CalModel.xaxislabels(Model.CalModel.selected_vars(1));
           b=Model.CalModel.xaxislabels(Model.CalModel.selected_vars(end));
       else
           a=Model.CalModel.xaxislabels(Model.CalModel.allint(Model.CalModel.selected_intervals,2));
           b=Model.CalModel.xaxislabels(Model.CalModel.allint(Model.CalModel.selected_intervals,3));
       end
       for i=1:max(size(a))
          h1temp=area([a(i) b(i)],[actualaxis(3) actualaxis(3)]); % Negative areas
          set(h1temp,'FaceColor',[0.75 0.75 0.75]);
          h2temp=area([a(i) b(i)],[actualaxis(4) actualaxis(4)]);
          set(h2temp,'FaceColor',[0.75 0.75 0.75]);
       end
		plot(Model.CalModel.xaxislabels,Model.CalModel.rawX,'k') % To overlay spectra on area plot
	end
	hold off
end % switch Model.type

figure('Position',pos2);
% Predicted versus measured for combined intervals
switch Model.type
  case 'iPLS'
	if strcmp(Model.val_method,'test')
        plotYref=Model.rawY(Model.segments,y_variable);
        plotYpred=Model.PLSmodel{interval}.Ypred(Model.segments,y_variable,no_of_int_lv);
        samplelabels=num2str(Model.segments);
    else
        plotYref=Model.rawY(:,y_variable);
        plotYpred=Model.PLSmodel{interval}.Ypred(:,y_variable,no_of_int_lv);
        samplelabels=num2str((1:size(plotYref,1))');
    end
    plot(plotYref,plotYpred,'w')
	text(plotYref,plotYpred,samplelabels);
	a=min([plotYref;plotYpred]);
	b=max([plotYref;plotYpred]);
	axis([a-abs(a)*0.1 b+abs(b)*0.1 a-abs(a)*0.1 b+abs(b)*0.1]);
	s=sprintf([titletext ', with %g PLS comp. for y-var. no. %g'],no_of_int_lv,y_variable);
	title(s);
	xlabel('Measured');
	ylabel('Predicted');
	hold on
	plot([a-abs(a)*0.1 b+abs(b)*0.1],[a-abs(a)*0.1 b+abs(b)*0.1]);
    %plot([a*0.9 b*1.1],[a*0.9 b*1.1]);
	hold off
	
	r=corrcoef([plotYref plotYpred]);
	s1 = sprintf('r = %0.4f',r(1,2));
	[RMSE,Bias]=rmbi(plotYref,plotYpred);

    if strcmp(lower(Model.val_method),'test')
        s2 = sprintf('RMSEP = %0.4f',RMSE);
    else
        s2 = sprintf('RMSECV = %0.4f',RMSE);        
    end
    s3 = sprintf('Bias = %0.4f',Bias);
    text(a+abs(a*0.08),1.1*b-abs(b*0.05),s1)
    text(a+abs(a*0.08),1.1*b-abs(b*0.10),s2)
    text(a+abs(a*0.08),1.1*b-abs(b*0.15),s3)
       
  case {'PLS'}
	if strcmp(Model.val_method,'test')
        plotYref=Model.rawY(Model.segments,y_variable);
        plotYpred=Model.PLSmodel{1}.Ypred(Model.segments,y_variable,no_of_int_lv);
       samplelabels=num2str(Model.segments);
    else
        plotYref=Model.rawY(:,y_variable);
        plotYpred=Model.PLSmodel{1}.Ypred(:,y_variable,no_of_int_lv);
        samplelabels=num2str((1:size(plotYref,1))');
    end
    plot(plotYref,plotYpred,'w')
	text(plotYref,plotYpred,samplelabels);
	a=min([plotYref;plotYpred]);
	b=max([plotYref;plotYpred]);
	axis([a-abs(a)*0.1 b+abs(b)*0.1 a-abs(a)*0.1 b+abs(b)*0.1]);
	s=sprintf([titletext ', with %g PLS comp. for y-var. no. %g'],no_of_int_lv,y_variable);
	title(s);
	xlabel('Measured');
	ylabel('Predicted');
	hold on
	plot([a-abs(a)*0.1 b+abs(b)*0.1],[a-abs(a)*0.1 b+abs(b)*0.1]);
	hold off
	
	r=corrcoef([plotYref plotYpred]);
	s1 = sprintf('r = %0.4f',r(1,2));
	[RMSE,Bias]=rmbi(plotYref,plotYpred);

    if strcmp(lower(Model.val_method),'test')
        s2 = sprintf('RMSEP = %0.4f',RMSE);
    else
        s2 = sprintf('RMSECV = %0.4f',RMSE);        
    end
	s3 = sprintf('Bias = %0.4f',Bias);
	text(a+abs(a*0.08),1.1*b-abs(b*0.05),s1)
    text(a+abs(a*0.08),1.1*b-abs(b*0.10),s2)
    text(a+abs(a*0.08),1.1*b-abs(b*0.15),s3)
    
  case {'PLSprediction'}
    plotYref=Model.Yref(:,y_variable);
    plotYpred=Model.Ypred(:,y_variable,no_of_int_lv);
   	plot(plotYref,plotYpred,'w')
	text(plotYref,plotYpred,num2str((1:size(plotYref,1))'));
	a=min([plotYref;plotYpred]);
	b=max([plotYref;plotYpred]);
	axis([a-abs(a)*0.1 b+abs(b)*0.1 a-abs(a)*0.1 b+abs(b)*0.1]);
	s=sprintf([titletext ', with %g PLS comp. for y-var. no. %g'],no_of_int_lv,y_variable);
	title(s);
	xlabel('Measured');
	ylabel('Predicted');
	hold on
	plot([a-abs(a)*0.1 b+abs(b)*0.1],[a-abs(a)*0.1 b+abs(b)*0.1]);
	hold off
	
	r=corrcoef([plotYref plotYpred]);
	s1 = sprintf('r = %0.4f',r(1,2));
	[RMSE,Bias]=rmbi(plotYref,plotYpred);
	s2 = sprintf('RMSEP = %0.4f',RMSE);
	s3 = sprintf('Bias = %0.4f',Bias);
	text(a+abs(a*0.08),1.1*b-abs(b*0.05),s1)
    text(a+abs(a*0.08),1.1*b-abs(b*0.10),s2)
    text(a+abs(a*0.08),1.1*b-abs(b*0.15),s3)
end

function [RMSE,Bias]=rmbi(Yref,Ypred)
[n,m]=size(Yref);
RMSE = sqrt( sum(sum((Ypred-Yref).^2))/(n*m) );
Bias = sum(sum(Ypred-Yref))/(n*m);
