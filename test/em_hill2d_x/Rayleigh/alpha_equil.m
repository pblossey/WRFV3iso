function [alpha] = alpha_equil(T,isotope,phase)

% function [alpha] = alpha_equil(T,isotope,phase)
% equilibrium fractionation for exchange between ice/liquid and
% vapor.  The arguments:
%
%  - T, absolute temperature in K
%  - isotope = 18 for H2O, 19 for HDO and 20 for H2O18.
%  - phase = 'liq' for liquid and 'ice' for ice.
%
%  e.g.  alpha_equil(280,19,'liq') will give the equilibrium
%  fractionation of HDO in vapor over liquid.

if not(strcmp(phase,'liq')) & not(strcmp(phase,'ice'))
  error(['Error: phase argument of alpha_equil can only have values ''liq'' ' ...
         'or ''ice''.'])
end

if isotope == 18
  % standard isotope, no fractionation
    coef = [0 0 0];
elseif isotope == 19
  % HDO
  if strcmp(phase,'liq')
    coef = [24884 -76.248 0.052612];
  else
    coef = [16288 0 -0.0934];
  end
elseif isotope == 20
  % H2O18
  if strcmp(phase,'liq')
    coef = [1137 -0.4156 -0.0020667];
  else
    coef = [0 11.839 -0.028224];
  end
else
  error('unrecognized isotope in alpha_equil');
end
  
alpha = exp( (coef(1) + T.*(coef(2) + T*coef(3))) ./ T.^2 );
