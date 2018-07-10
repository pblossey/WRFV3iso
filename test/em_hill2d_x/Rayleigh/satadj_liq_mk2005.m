function [T] = satadj_liq_mk2005(h,qt,z,p)

% FUNCTION [T] = SATADJ(H,QT,Z,P) computes the temperature of moist air
% that is induced by the condensation/evaporation of liquid water (NO
% ICE PERMITTED), resulting in a mixture that is at equilibrium.  The
% equilibrium mixture will either have no liquid water or its water
% vapor mixing ratio will be equal to the saturation value.

global Cp L g

T0 = (h - g*z - L*qt)/Cp;
if qsw(p,T0)>=qt
  T=T0;
else
  T = fzero(@(T) testsat(T,h,qt,z,p),T0);
end

function [residual] = testsat(T,h,qt,z,p)

global Cp L g

residual = h - (Cp*T + g*z + L*min(qt,qsw_mk2005(p,T)));

