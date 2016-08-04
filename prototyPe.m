
load inputdata;
interim1 = cov(y_train);
interim2 = cov(y_train) * banding(size(y_train,1),1);

y_train_mod = inv(sqrt(interim2)) * y_train;
x_train_mod = inv(sqrt(interim2)) * x_train;

load inputdata;
[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,10,'CV',10);
yfit_PLS_train = [ones(size(x_train,1),1) x_train]*betaPLS;

[XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train_mod,10,'CV',10);
yfit_PLS_train_mod = [ones(size(x_train,1),1) x_train]*betaPLS;

corr(yfit_PLS_train,y_train)
corr(yfit_PLS_train_mod,y_train_mod)