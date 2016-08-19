%collate results

if exist('D:\NIR Gui Project\results_varSelection\uvepls_record.mat','file') ~= 2
return
end

if exist('D:\NIR Gui Project\var_file.xlsx','file') == 2
delete('D:\NIR Gui Project\var_file.xlsx');
end

load uvepls_record.mat
filename = 'var_file.xlsx';
Method = {'MC-UVE-PLS'};

T1 = table(Method,C_uvepls, R2Train_uvepls, RMSECV_uvepls, R2Test_uvepls, RMSEP_uvepls, Error_Percent_uvepls);
%renaming the headers for clarity
T1.Properties.VariableNames {'C_uvepls'} = 'C';
T1.Properties.VariableNames {'R2Train_uvepls'} = 'R2Train';
T1.Properties.VariableNames {'RMSECV_uvepls'} = 'RMSECV';
T1.Properties.VariableNames {'R2Test_uvepls'} = 'R2Test';
T1.Properties.VariableNames {'RMSEP_uvepls'} = 'RMSEP';
T1.Properties.VariableNames {'Error_Percent_uvepls'} = 'Error_Percent';

writetable(T1,filename,'Sheet',1,'Range','A1');

if exist('D:\NIR Gui Project\results_varSelection\bipls_record.mat','file') == 2
    load bipls_record.mat
    T2 = [{'biPLS'},C_bipls,R2Train_bipls, RMSECV_bipls, R2Test_bipls, RMSEP_bipls, Error_Percent_bipls];
    xlswrite(filename,T2,1,'A3');
end

if exist('D:\NIR Gui Project\results_varSelection\cars_record.mat','file') == 2
    load cars_record.mat
    T3 = [{'CARS'},C_cars, R2Train_cars, RMSECV_cars, R2Test_cars, RMSEP_cars, Error_Percent_cars];
    xlswrite(filename,T3,1,'A4');
end

if exist('D:\NIR Gui Project\results_varSelection\gapls_record.mat','file') == 2
    load gapls_record.mat
    T4 = [{'GA-PLS'},C_gapls, R2Train_gapls, RMSECV_gapls, R2Test_gapls, RMSEP_gapls, Error_Percent_gapls];
    xlswrite(filename,T4,1,'A5');
end

if exist('D:\NIR Gui Project\results_varSelection\jkpls_record.mat','file') == 2
    load jkpls_record.mat
    T5 = [{'JK-PLS'},C_jkpls, R2Train_jkpls, RMSECV_jkpls, R2Test_jkpls, RMSEP_jkpls, Error_Percent_jkpls];
    xlswrite(filename,T5,1,'A6');
end

if exist('D:\NIR Gui Project\results_varSelection\vcppls_record.mat','file') == 2
    load vcppls_record.mat
    T6 = [{'VCPA-PLS'},C_vcppls, R2Train_vcppls, RMSECV_vcppls, R2Test_vcppls, RMSEP_vcppls, Error_Percent_vcppls];
    xlswrite(filename,T6,1,'A7');
end

fclose('all');

winopen('var_file.xlsx')
clear all;
