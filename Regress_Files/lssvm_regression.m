function [ r2train,r2test,rmsecvC,rmsepC  ] = lssvm_regression( x_train, y_train, x_test , y_test )
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

addpath LSSVM_files

%% Model Construction

% Construct test data
Xt = x_test;
Yt = y_test;

% training data
X = x_train;
Y = y_train;

h = waitbar(0);
model = initlssvm(X,Y,'f',[],[],'RBF_kernel','o')
waitbar(1 / 3);
model = tunelssvm(model,'simplex','crossvalidatelssvm',{10,'mse'})
waitbar(2 /3);
model = trainlssvm(model)
close(h)
%Predict the training values
Yhci = simlssvm(model,X); % Add K-cross CV and get tuned results or so !
%Predict the test values
Yhpi = simlssvm(model,Xt);

corr_train = corr(Yhci,y_train)
figure
plot(y_train,Yhci,'o')
corr_test = corr(Yhpi,y_test)
figure
plot(y_test,Yhpi,'o')

s1 = sum((y_train - Yhci).^2);
s2 = sum((y_train - mean(y_train)).^2);
r2train = (1 - (s1/s2))*100
rmsecvC = sqrt(mean((y_train-Yhci).^2))

s11 = sum((y_test - Yhpi).^2);
s22 = sum((y_test - mean(y_test)).^2);
r2test = (1 - (s11/s22))*100
rmsepC = sqrt(mean((y_test-Yhpi).^2))

Error_Percent = abs(mean(abs(Yhpi-y_test)))/((max(y_train)-min(y_train)))*100

% Printing the values to respective matfile
%C = find(r2test == max(r2test)) ;
R2Train_lssvm = r2train;
R2Test_lssvm =r2test;
RMSECV_lssvm =rmsecvC;
RMSEP_lssvm =rmsepC;
Error_Percent_lssvm =Error_Percent;

if exist('D:\NIR Gui Project\results_regSelection\lssvm_record.mat','file') == 2
delete('D:\NIR Gui Project\results_regSelection\lssvm_record.mat');
else
save('D:\NIR Gui Project\results_regSelection\lssvm_record.mat','R2Train_lssvm','RMSECV_lssvm','R2Test_lssvm','RMSEP_lssvm','Error_Percent_lssvm');
end


end

