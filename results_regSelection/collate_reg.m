%collate results

if exist('D:\NIR Gui Project\results_regSelection\pls_record.mat','file') ~= 2
return
end

if exist('D:\NIR Gui Project\reg_file.xlsx','file') == 2
delete('D:\NIR Gui Project\reg_file.xlsx');
end

load pls_record.mat
filename = 'reg_file.xlsx';
Method = {'pls'} ;

T1 = table(Method,C_pls, R2Train_pls, RMSECV_pls, R2Test_pls, RMSEP_pls, Error_Percent_pls);
%renaming the headers for clarity
T1.Properties.VariableNames {'C_pls'} = 'C';
T1.Properties.VariableNames {'R2Train_pls'} = 'R2Train';
T1.Properties.VariableNames {'RMSECV_pls'} = 'RMSECV';
T1.Properties.VariableNames {'R2Test_pls'} = 'R2Test';
T1.Properties.VariableNames {'RMSEP_pls'} = 'RMSEP';
T1.Properties.VariableNames {'Error_Percent_pls'} = 'Error_Percent';

writetable(T1,filename,'Sheet',1,'Range','A1');

if exist('D:\NIR Gui Project\results_regSelection\lssvm_record.mat','file') == 2
    load lssvm_record.mat
    T2 = [{'lssvm'},0,R2Train_lssvm, RMSECV_lssvm, R2Test_lssvm, RMSEP_lssvm, Error_Percent_lssvm];
    xlswrite(filename,T2,1,'A3');
end

if exist('D:\NIR Gui Project\results_regSelection\pca_record.mat','file') == 2
    load pca_record.mat
    T3 = [{'pcr'},C_pca, R2Train_pca, RMSECV_pca, R2Test_pca, RMSEP_pca, Error_Percent_pca];
    xlswrite(filename,T3,1,'A4');
end

if exist('D:\NIR Gui Project\results_regSelection\ridge_record.mat','file') == 2
    load ridge_record.mat
    T4 = [{'ridge'},C_ridge, R2Train_ridge, RMSECV_ridge, R2Test_ridge, RMSEP_ridge, Error_Percent_ridge];
    xlswrite(filename,T4,1,'A5');
end

if exist('D:\NIR Gui Project\results_regSelection\lasso_record.mat','file') == 2
    load lasso_record.mat
    T5 = [{'lasso'},C_lasso, R2Train_lasso, RMSECV_lasso, R2Test_lasso, RMSEP_lasso, Error_Percent_lasso];
    xlswrite(filename,T5,1,'A6');
end

fclose('all');

winopen('reg_file.xlsx')
clear all;
