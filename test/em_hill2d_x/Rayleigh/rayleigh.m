clear t theta q p R logR

t0 = 300;
theta0 = 300;
q0 = 20e-3;
p0 = 1e3;
R0 = 0.9;

t(1) = t0;
theta(1) = theta0;
q(1) = q0;
p(1) = p0;
R(1) = R0;
logR(1) = log(R0);

dp = 1;
for n = 1:950/dp
  p(n+1) = p(n) - dp;
  theta(n+1) = theta(n);
  t(n+1) = theta(n+1)*(p(n+1)/1000).^(Rd/Cp);
  q(n+1) = q(n);

  if t(n+1)<0
    tnew = satadj_ice_mk2005(Cp*t(n+1) + L*q(n+1),q(n+1),0,100*p(n+1));
    t(n+1) = tnew;
    q(n+1) = min(q(n+1),qsi_mk2005(100*p(n+1),t(n+1)));
    logR(n+1) = logR(n) + (alpha_equil_ice(t(n+1),19)-1) ...
        *(log((28/18)*q(n+1)) - log((28/18)*q(n)));
  else
    tnew = satadj_liq_mk2005(Cp*t(n+1) + L*q(n+1),q(n+1),0,100*p(n+1));
    t(n+1) = tnew;
    q(n+1) = min(q(n+1),qsw_mk2005(100*p(n+1),t(n+1)));
    logR(n+1) = logR(n) + (alpha_equil_liq(t(n+1),19)-1) ...
        *(log((28/18)*q(n+1)) - log((28/18)*q(n)));
  end

  theta(n+1) = t(n+1)*(1000/p(n+1)).^(Rd/Cp);

  R(n+1) = exp(logR(n+1));

  if n==1 | mod(p(n+1),50)==0
    disp([p(n+1) t(n+1) theta(n+1) 1e3*q(n+1) R(n+1)])
  end
end
