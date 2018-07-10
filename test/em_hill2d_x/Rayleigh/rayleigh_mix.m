function [qvmix,dDmix,dOmix] = rayleigh_mix(qv,dD,dO,dD1,dD2)

%  function [qvmix,dDmix,dOmix] = rayleigh_mix(qv,dD,dO,dD1,dD2)
%
% compute mixing curves of qv, deltaD and deltaO18 between two
% parcels on a given Rayleigh curve (or any fractionation curve)
% identified by two deltaD values.

    rind = find(dD<max(dD));
                                                     
    qv1 = interp1(dD(rind),qv(rind),dD1);
    dO1 = interp1(dD(rind),dO(rind),dD1);
    qD1 = (1 + dD1/1000)*qv1;
    qO1 = (1 + dO1/1000)*qv1;

    qv2 = interp1(dD(rind),qv(rind),dD2);
    dO2 = interp1(dD(rind),dO(rind),dD2);
    qD2 = (1 + dD2/1000)*qv2;
    qO2 = (1 + dO2/1000)*qv2;

    tmp = linspace(0,1);
    qvmix = qv1*tmp + qv2*(1-tmp);
    qDmix = qD1*tmp + qD2*(1-tmp);
    qOmix = qO1*tmp + qO2*(1-tmp);
    dDmix = 1000*(qDmix./qvmix - 1);
    dOmix = 1000*(qOmix./qvmix - 1);
                                                   
