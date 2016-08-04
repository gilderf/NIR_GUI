function [ Y ] = property_selector( Y_complex,prop_selected )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch prop_selected{1}
    case 'Cart. Thickness'
        Y = Y_complex.thickness;
    case 'Cart. Cal. Thickness'
        Y = Y_complex.thickness_cal;
    case 'Inst.Mod'
        Y = Y_complex.instant;
    case 'Dynamic.Mod'
        Y = Y_complex.freq1;
    case 'Equilibrium.Mod'
        Y = Y_complex.equilibrium;
    otherwise
        warning('Invalid Property')
end

% Y_complex.thickness = myRef(:,1);
% Y_complex.thickness_cal = myRef(:,2);
% Y_complex.ICRS = myRef(:,3);
% Y_complex.instant = myRef(:,5);
% Y_complex.instant = Y_complex.instant./1e6;
% Y_complex.equilibrium = myRef(:,7);
% Y_complex.equilibrium = Y_complex.equilibrium./1e6;
% Y_complex.freq1 = myRef(:,8);
% Y_complex.freq1 = freq1./1e6;

end

