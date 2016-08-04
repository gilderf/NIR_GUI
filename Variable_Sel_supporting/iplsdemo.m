echo off
close all
clc
disp(' ')
disp('  This is a short demo of iPLS in the iToolbox')
disp(' ')
disp('  We want to make an iPLS model of near infrared data measured')
disp('  on 40 beers with original extract as y-variable')
disp('  (20 samples are in an independent test set not used in the demo)')
disp(' ')
disp('  Load the data')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
load nirbeer; echo off
disp(' ')
disp('  The model is run with 20 intervals, 10 PLS components, mean centered data,')
disp('  and 5 segments using syst123 cross validation')
disp(' ')
disp('  Use Model=ipls(X,Y,no_of_lv,prepro_method,intervals,xaxislabels,val_method,segments)')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
Model=ipls(Xcal,ycal,10,'mean',20,xaxis,'syst123',5); echo off
disp(' ')
disp('  Make the iPLS plot:')
disp('  Use iplsplot(Model,labeltype,optimal_lv_global,max_yaxis,plottype)')
disp('  optimal_lv_global,max_yaxis,plottype are optional')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
iplsplot(Model,'intlabel'); echo off
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
close all

disp('  Make iPLS plot with wavelengths as labels')
disp('  Press any key to continue')
disp(' ')
pause
echo on
iplsplot(Model,'wavlabel'); echo off
disp('  Press any key to continue')
disp(' ')
pause
close all

disp('  Interval number 10 has the best performance in this case')
disp(' ')
disp('  Plot predicted versus measured for interval 10')
disp('  with 6 PLS components')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
plspvsm(Model,6,10); echo off
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
close all

disp('  Plot RMSECV for the global model')
disp('  Use plsrmse(Model,interval) with interval=0')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
plsrmse(Model,0); echo off
disp(' ')
disp('  Due to the noisy part in the spectra a minimum is found at 1 PLS components') 
disp('  but it is more fair to compare the interval model with e.g. a 4 PLS components model') 
disp('  Press any key to continue')
pause
disp(' ')
echo on
iplsplot(Model,'intlabel',4); echo off
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
close all

disp('  Plot predicted versus measured plot for the global model with 4 PLS components')
disp('  for comparison')
disp('  Use plspvsm(Model,no_of_int_lv,interval,y_variable) with interval=0')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
plspvsm(Model,4,0,1); echo off
disp(' ')
disp('  Press any key to continue')
pause
disp(' ')
close all

disp('  If you want to see the actual intervals and/or the corresponding')
disp('  wavelength variables use intervals(Model)')
disp(' ')
echo on
intervals(Model); echo off
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
disp(' ')
disp('  All model information is stored in Model, which is a structure array. Just type')
disp('  Model at the command line to display the contents')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
Model
echo off
disp(' ')
disp('  Also try Model.PLSmodels{10} to see the outputs from the PLS model on the 10th interval')
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
echo on
Model.PLSmodel{10}
echo off
disp(' ')
disp('  Press any key to continue')
disp(' ')
pause
disp('  You can now proceed with ordinary test set prediction (plspredict.m)')
disp('  to see if the selected interval gives comparable results for an independent test set')
disp(' ')
disp('  END OF DEMO')
