function [alpha] = alpha_equil_ice(T,isotope)

% function [alpha] = alpha_equil_ice(T)
% equilibrium fractionation for exchange between ice and vapor.

if isotope == 18
  % standard isotope, no fractionation
  coef = [0 0 0];
elseif isotope == 19
  % HDO
  coef = [16288 0 -0.0934];
elseif isotope == 20
  % H2O18
  coef = [0 11.839 -0.028224];
else
  error('unrecognized isotope');
end
  
alpha = exp( (coef(1) + T.*(coef(2) + T*coef(3))) ./ T.^2 );
