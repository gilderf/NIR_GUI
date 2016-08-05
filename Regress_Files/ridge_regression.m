function [ r2train,r2test,rmsecvC,rmsepC  ] = ridge_regression( x_train, y_train, x_test , y_test , func_param)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

interim = [0 : .1 : 1000];
%interim = [0 : .001 : 1];
till = size(interim,2);

r2train = zeros(till,1);
rmsecvC = zeros(till,1);
r2test = zeros(till,1);
rmsepC = zeros(till,1);

Error_Percent = zeros(till,1);

%C = zeros(till,1);
h = waitbar(0);
for i = 1:till
    waitbar(i / till);
    
    [B] = ridge(y_train,x_train,interim(i),0);
    
    y_predicted_train = B(1) + x_train*B(2:end);
    y_predicted_test = B(1) + x_test*B(2:end);
    
    s1 = sum((y_train - y_predicted_train).^2);
    s2 = sum((y_train - mean(y_predicted_train)).^2);
    r2train(i) = (1 - (s1/s2))*100;
    rmsecvC(i) = sqrt(mean((y_train-y_predicted_train).^2));
    
    s11 = sum((y_test - y_predicted_test).^2);
    s22 = sum((y_test - mean(y_predicted_test)).^2);
    r2test(i) = (1 - (s11/s22))*100;
    rmsepC(i) = sqrt(mean((y_test-y_predicted_test).^2));
    
    Error_Percent(i) = abs(mean(abs(y_predicted_test-y_test)))/(max(y_train)-min(y_train))*100;
    
    
end
close(h)

% Print the results Grid
ncompmax = till;
figure;
title('Baseline per manuscript');
subplot(2,3,1);
plot(1:ncompmax,r2train)
xlabel('Lambda Value');
ylabel('R^2 Train');
hold on
plot(find(r2train == min(r2train)),min(r2train),'o')
plot(find(r2test == max(r2test)),r2train,'X')
hold off

subplot(2,3,2);
plot(1:ncompmax,r2test)
xlabel('Lambda Value');
ylabel('R^2 Test');
hold on
plot(find(r2test == max(r2test)),max(r2test),'o')
hold off

subplot(2,3,3);
plot(1:ncompmax,rmsepC)
xlabel('Lambda Value');
ylabel('rmsepC');
hold on
plot(find(rmsepC == min(rmsepC)),min(rmsepC),'o')
plot(find(r2test == max(r2test)),rmsepC,'X')
hold off

subplot(2,3,4);
plot(1:ncompmax,rmsecvC)
xlabel('Lambda Value');
ylabel('RMSECV');
hold on
plot(find(rmsecvC == min(rmsecvC)),min(rmsecvC),'o')
plot(find(r2test == max(r2test)),rmsecvC,'X')
hold off
title_to_display = strcat('Ridge Regression : ',func_param);
suptitle(title_to_display)

subplot(2,3,5);
plot(1:ncompmax,Error_Percent)
xlabel('Lambda Value');
ylabel('Error %');
hold on
plot(find(Error_Percent == min(Error_Percent)),min(Error_Percent),'o')
plot(find(r2test == max(r2test)),Error_Percent,'X')
hold off

% Printing the values to respective matfile
C_ridge = find(r2test == max(r2test)) ;
C = C_ridge;
R2Train_ridge = r2train(C);
R2Test_ridge =r2test(C);
RMSECV_ridge =rmsecvC(C);
RMSEP_ridge =rmsepC(C);
Error_Percent_ridge =Error_Percent(C);

%if there is a previous copy of the file delete it and make a new one
if exist('D:\NIR Gui Project\results_regSelection\ridge_record.mat','file') == 2
delete('D:\NIR Gui Project\results_regSelection\ridge_record.mat');
end
save('D:\NIR Gui Project\results_regSelection\ridge_record.mat','C_ridge','R2Train_ridge','RMSECV_ridge','R2Test_ridge','RMSEP_ridge','Error_Percent_ridge');

end

