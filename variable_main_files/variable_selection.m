function [ wavelength_final ] = variable_selection( x_train,y_train,x_test,y_test,var_sel_mthd,range )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

var_sel_mthd = var_sel_mthd{1};

switch var_sel_mthd
    case 'UVE'
        [ wavelength_final  ] = UVE_var_processor( x_train, y_train, x_test ,y_test);
    case 'CARS'
        [ wavelength_final  ] = CARS_var_processor( x_train, y_train, x_test ,y_test);
    case 'VCP'
        [ wavelength_final  ] = VCP_var_processor( x_train, y_train, x_test , y_test);
    case 'biPLS'
        [ wavelength_final  ] = biPLS_var_processor( x_train, y_train, x_test ,y_test, range);
    case 'GA'
        [ wavelength_final  ] = GA_var_processor( x_train, y_train, x_test ,y_test);
    case 'JK'
        [ wavelength_final  ] = Jack_Knife_var_processor( x_train, y_train,x_test,y_test);
    case 'RunAll'
        run_all_var_sel( x_train,y_train,x_test,y_test,var_sel_mthd,range);
    otherwise
        warning('no Regression selected')
        
        
end

