function [alpha] = alpha_equil_liq(T,isotope)

% function [alpha] = alpha_equil_liq(T)
% equilibrium fractionation for exchange between liquid and vapor.

if isotope == 18
  % standard isotope, no fractionation
  coef = [0 0 0];
elseif isotope == 19
  % HDO
  coef = [24884 -76.248 0.052612];
elseif isotope == 20
  % H2O18
  coef = [1137 -0.4156 -0.0020667];
else
  error('unrecognized isotope');
end
  
alpha = exp( (coef(1) + T.*(coef(2) + T*coef(3))) ./ T.^2 );
