function [T] = satadj(h,qt,z,p,phase,formula)

% FUNCTION [T] = SATADJ(H,QT,Z,P) computes the temperature of moist air
% that is induced by the condensation/evaporation of liquid water (NO
% ICE PERMITTED), resulting in a mixture that is at equilibrium.  The
% equilibrium mixture will either have no liquid water or its water
% vapor mixing ratio will be equal to the saturation value.

global Cp L g

if nargin < 6; formula = 'flatau'; end
if nargin < 5; phase = 'liq'; end

T0 = (h - g*z - L*qt)/Cp;
if qsat(p,T0,phase,formula)>=qt
  T=T0;
else
  T = fzero(@(T) testsat(T,h,qt,z,p,phase,formula),T0);
end

function [residual] = testsat(T,h,qt,z,p,phase,formula)

global Cp L g

residual = h - (Cp*T + g*z + L*min(qt,qsat(p,T,phase,formula)));

