function Result=VCPA(X,Y,A_max,fold,method,BMS_Run,EDF_Run,Ratio_Better)

%+++ Variable Combination population Analysis for Variable selection.
%+++ X: The data matrix of size m x n
%+++ y: The reponse vector of size m x 1
%+++ A_max: the maximal principle to extract.
%+++ fold: the group number for cross validation.
%+++ method: pretreatment method.
%+++ BMS_Run: the number of binary matrix sampling(BMS) run
%+++ EDF_Run: the number of exponentially decreasing function(EDF) run
%+++ Ratio_Better: Ratio of better models from a population model 
%+++ Yonghuan Yun, Dec.15, 2013. yunyonghuan@foxmail.com
%+++ Advisor: Yizeng Liang, yizeng_liang@263.net
%+++ If you want Latest version of VCPA, please contact
%+++ yunyonghuan@foxmail.com.
%+++ Reference: Using Variable Combination Population Analysis for 
%               Variable Selection in Multivariate Calibration£¬Analytica Chimica Acta£¬
%               doi:10.1016/j.aca.2014.12.048

% This version of VCPA inclines to select fewer variables because it
% finds the optimal variable subset from all combinations of the remaining L variables
% after many EDF runs. L can be set to larger value which depends on your computer.

if nargin<8;Ratio_Better=0.1;end;
if nargin<7;EDF_Run=50;end;
if nargin<6;BMS_Run=1000;end;
if nargin<5;method='center';end;
if nargin<4;fold=5;end;
if nargin<3;A_max=10;end;

XX=X;
[~,Nx]=size(XX);
[~,N]=size(X);   % N will change in EDF loop
VarIndex=1:N;       
VarID=1:Nx;
Ratio=zeros(1,EDF_Run);

% This parameters can be set to larger, which depends on your computer.
L=14;% the number of left variables for the final run. It can also be 15, 16, 17, 18... which depends on your computer.

%+++ exponentially decreasing function (EDF)
r0=1;  r1=L/Nx;  b=log(r0/r1)/(EDF_Run);  a=1;             %+++ Parameter of exponentially decreasing function. 
for iter=1:EDF_Run
     
     ratio=a*exp(-b*(iter));                      %+++ Ratio of retained variables through emliminating variables by force.
     Ratio(iter)=ratio;
       
end

%+++ the mean number of variables in every BMS sampling. This step can be also tuned by user. 
% sqrt(N) can be changed as other form which depneds on users. 
predefinevar=round(linspace(floor(sqrt(N)),L/2,EDF_Run)); 
frequency=zeros(iter,Nx);

%+++ Main Loop
tic;
h = waitbar(0);
for iter=1:EDF_Run
    waitbar(iter / EDF_Run);
     binary_matrix=generate_binary_matrix(N, BMS_Run,predefinevar(iter));  %+++ generate binary matrix for sampling in the varibale space
     for i=1:BMS_Run  
         temp = binary_matrix(i,:); 
         del_X = temp==0;  
         X_new = X;   
         X_new(:,del_X) = [];
         CV=plscvfold(X_new,Y,A_max,fold,method);   
         RMSECV(i)=CV.RMSECV;   
     end    
     
     VarNum_retain=round(Nx*Ratio(iter));            %+++ the number of retained variables in each EDF run
     [~,rank_CV]=sort(RMSECV);      %+++ Sort the RMSECV  
     best_matrix=binary_matrix(rank_CV(1:BMS_Run*Ratio_Better),:);%+++ Choose the better models
     [~,rank_variables]=sort(sum(best_matrix),'descend');        %+++ Sort the variable based on their frequency in the better models
     f=zeros(1,Nx);
     f(VarIndex)=sum(best_matrix);
     frequency(iter,:)=f;   
     VarIndex=sort(VarIndex(rank_variables(1:VarNum_retain)));              %+++ Retain the variable based on EDF
     X=XX(:,VarIndex);
     [~,N]=size(X); 
     minRMSECV_EDF(iter)=min(RMSECV);
     clear binary_matrix RMSECV best_matrix
     if rem(iter,10)==0;  
         fprintf('The %d(th)/%d EDF running has finished,elapsed time is %g seconds!!\n',iter,EDF_Run,toc);    %+++ Screen output.
     end
end
close(h)

%+++ Gennerate all combination of the left varibales
l=2^length(VarIndex)-1;
All_Combintaion=zeros(l,Nx);
temp_zeros=zeros(1,Nx);
k=1;
for i=1:length(VarIndex)
    AllCom_temp{i}=combntns(VarIndex,i);            
    for j=1:size(AllCom_temp{i},1)        
        temp_zeros(AllCom_temp{1,i}(j,:))=1;
        All_Combintaion(k,:)=temp_zeros;
        temp_zeros=zeros(1,Nx);
        k=k+1;
    end

end


%+++ The final run for computing the RMSECV of all combination of left L varibales
for i=1:size(All_Combintaion,1) 
         temp = All_Combintaion(i,:); 
         del_X = temp==0;  
         X_new = XX;   
         X_new(:,del_X) = [];
         CV=plscvfold(X_new,Y,A_max,fold,method);   
         RMSECV_final(i)=CV.RMSECV;   
            if rem(i,1000)==0;        
                fprintf('The %d(th)/%d of the final run has finished,elapsed time is %g seconds!!\n', i,size(All_Combintaion,1),toc);       
            end
end
    
[~,Rank_final]=min(RMSECV_final);

selected_variables=VarID(All_Combintaion(Rank_final,:)==1);  %+++ select the variables of the lowest RMSECV in the final run
CV=plscvfold(XX(:,selected_variables),Y,A_max,fold,method);   


 Result.A_max=A_max;
 Result.fold=fold;
 Result.Ratio=Ratio; %+++ The ratio of each EDF run
 Result.minRMSECV=min(CV.RMSECV);%+++ The RMSECV of the selected variables
 Result.minRMSECV_EDF=minRMSECV_EDF;%+++ The minimum RMSECV of each EDF run
 Result.Q2_Vsel=CV.Q2_max;%+++ The Q2 of the selected variables
 Result.Vsel=selected_variables;




