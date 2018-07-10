function [p,theta,qv,delHDO,delO18] = rayleigh_curve(theta0,p0,q0,delHDO0,delO180,T00)

% function [p,theta,qv,delHDO,delO18] = rayleigh_curve(theta0,p0,q0,delHDO0,delO180,T00)
%
%  computes the Rayleigh curve assuming the surface parcel conditions
%  specified in theta0 (K), p0 (Pa), q0 (kg/kg), delHDO (per mil),
%  delO18 (per mil) and, optionally, T00 (K).  The parcel is raised
%  adiabatically, and heat is released as vapor condenses.  The
%  optional argument T00 (default = 273.15 K), specifies the
%  temperature at which the saturation vapor pressure computation
%  changes from liquid to ice.  The condensate is removed from the
%  parcel immediately.  The isotopes HDO and O18 in the remaining
%  vapor fractionate as:
%
%    d(ln R) = (alpha - 1) d(ln q)
%
%  where R=[isotope]/[H2O]

global Rd Cp L

if nargin < 6
  T00 = 273.15;
end

t(1) = theta0*(p0/1e5)^(Rd/Cp);
theta(1) = theta0;
q(1) = q0;
p(1) = p0;
R(1,:) = 1+[delHDO0 delO180]/1000
logR(1,:) = log(R(1,:));

p = linspace(p0,1,round((p0-1)/100));
disp(p0)
disp(p(1:10))
dp = diff(p(1:2));
for n = 1:length(p)-1
% $$$   p(n+1) = p(n) - dp;
  theta(n+1) = theta(n);
  t(n+1) = theta(n+1)*(p(n+1)/1e5).^(Rd/Cp);
  q(n+1) = q(n);

  if t(n+1)<T00
    tnew = satadj(Cp*t(n+1) + L*q(n+1),q(n+1),0,p(n+1),'ice','mk2005');
    t(n+1) = tnew;
    q(n+1) = min([q(n+1) qsat(p(n+1),t(n+1),'ice','mk2005')]);
    logR(n+1,:) = logR(n,:) ...
        + ([alpha_equil(t(n+1),19,'ice') alpha_equil(t(n+1),20,'ice')]-1) ...
        *(log(q(n+1)) - log(q(n)));
  else
    tnew = satadj(Cp*t(n+1) + L*q(n+1),q(n+1),0,p(n+1),'liq','mk2005');
    t(n+1) = tnew;
    q(n+1) = min([q(n+1) qsat(p(n+1),t(n+1),'liq','mk2005')]);
    logR(n+1,:) = logR(n,:) ...
        + ([alpha_equil(t(n+1),19,'liq') alpha_equil(t(n+1),20,'liq')]-1) ...
        *(log(q(n+1)) - log(q(n)));
  end

  theta(n+1) = t(n+1)*(1e5/p(n+1)).^(Rd/Cp);
  R(n+1,:) = exp(logR(n+1,:));

  if n==1 | mod(p(n+1),3e3)<dp/2
    disp([p(n+1) t(n+1) theta(n+1) 1e3*q(n+1) R(n+1,:)])
  end
end

qv = q;
delHDO = 1000*(exp(logR(:,1))-1);
delO18 = 1000*(exp(logR(:,2))-1);
