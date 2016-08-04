%+++ Give an example of using VCPA  

load('Corn starch.mat')



A_max=10;
fold=5;
method='center';
Ratio_Better=0.1;
EDF_Run=50;
BMS_Run=1000;


Result=VCPA(Xtrain,Ytrain,A_max,fold,method,BMS_Run,EDF_Run,Ratio_Better);
selected_variables=Result.Vsel;
F=predict(Xtrain,Ytrain,Xtest,Ytest,selected_variables,A_max,fold,method);

% plot of EDF

plot(1:EDF_Run,Result.Ratio)

% plot of RMSECV of EDF run

plot(1:EDF_Run,Result.minRMSECV_EDF)

% hist the distribution of  mean number of the selected variables by BMS

predefinevar=round(sqrt(size(Xtrain,2)));
binary_matrix=generate_binary_matrix(N, BMS_Run,predefinevar);
bar(sum(binary_matrix)) % every variable has the same chance to be sampled 
hist(sum(binary_matrix'))


% 50 replicate running
tic;
for i=1:50
Result=VCPA(Xtrain,Ytrain,A_max,fold,method,BMS_Run,EDF_Run,Ratio_Better);
selected_variables=Result.Vsel;
F=predict(Xtrain,Ytrain,Xtest,Ytest,selected_variables,A_max,fold,method);
RMSEP(i)=F.RMSEP;
RMSEC(i)=F.RMSEC;
number(i)=length(selected_variables);
variables{i}=selected_variables;
RMSECV_EDF(i,:)=Result.minRMSECV_EDF;
fprintf('The %d(th) run has finished,elapsed time is %g seconds!!\n', i,toc);       
end

% plot of RMSECV of EDF run
plot(1:EDF_Run,mean(RMSECV_EDF))

