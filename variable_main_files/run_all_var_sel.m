function run_all_var_sel( x_train,y_train,x_test,y_test,var_sel_mthd,range )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% h = waitbar(0);
% waitbar(1 / 7);
UVE_var_processor( x_train, y_train, x_test ,y_test);
% waitbar(2 / 7);
CARS_var_processor( x_train, y_train, x_test ,y_test);
% waitbar(3 / 7);
biPLS_var_processor( x_train, y_train, x_test ,y_test, range);
% waitbar(4 / 7);
VCP_var_processor( x_train, y_train, x_test , y_test);
% waitbar(5 / 7);
Jack_Knife_var_processor( x_train, y_train,x_test,y_test);
% waitbar(6 / 7);
GA_var_processor( x_train, y_train, x_test ,y_test);
% waitbar(7 / 7);
% close(h)

end

