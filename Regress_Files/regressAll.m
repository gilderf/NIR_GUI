function regressAll( x_train, y_train, x_test, y_test, range,var_sel_mthd )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% h = waitbar(0);
% waitbar(1 / 7);
pls_regression( x_train, y_train, x_test , y_test, var_sel_mthd );
% waitbar(2 / 7);
pca_regression( x_train, y_train, x_test , y_test , range, var_sel_mthd);
% waitbar(3 / 7);
%ica_regression( x_train, y_train, x_test , y_test ,range);
% waitbar(4 / 7);
ridge_regression( x_train, y_train, x_test , y_test ,var_sel_mthd);
% waitbar(5 / 7);
lssvm_regression( x_train, y_train, x_test , y_test);
% waitbar(6 / 7);
lasso_regression( x_train, y_train, x_test , y_test, var_sel_mthd );
% waitbar(7 / 7);
% close(h)

end