clear all; define_constants

%%addpath('/home/disk/eos1/bloss/matlab/Rayleigh')
addpath('/n/home07/pblossey/matlab/Rayleigh')

in(1).nc = 'wrfinput_d01';

whoption = input(['Which option [1=No Fractionation, ' ...
                  '2=Rayliegh]: ']);

% copy the clean (isotope-free) WRF input file to save a clean copy
eval(sprintf('!cp %s %s_NoWISO',in(1).nc,in(1).nc))

switch whoption
 case 1
  % no fractionation
  wh = {'VAPOR','CLOUD','ICE','RAIN','SNOW','GRAUP'};
  whiso = {'HDO','O18'};
  for m = 1:length(wh)
    varname = sprintf('Q%s',wh{m})
    qq = ncread(in(1).nc,varname);
    for nn = 1:length(whiso)
      isoname = sprintf('%s_Q%s',whiso{nn},wh{m})
      ncwrite(in(1).nc,isoname,qq);
    end
  end

 case 2
  qv = ncread(in(1).nc,'QVAPOR');
  theta = 300+ncread(in(1).nc,'T');
  p = ncread(in(1).nc,'PB') + ncread(in(1).nc,'P'); 
  T = theta.*(p/1e5).^(Rd/Cp);

  rh = qv./qsw(p,T);

  figure(1); clf
  subplot(321); pcolor(squeeze(qv(:,1,:))'); shading flat; colorbar
  title('qv [kg kg^{-1}]');
  subplot(322); pcolor(squeeze(T(:,1,:))'); shading flat; colorbar
  title('T [K]');
  subplot(323); pcolor(squeeze(p(:,1,:))'); shading flat; colorbar
  title('p [Pa]');
  subplot(324); pcolor(squeeze(theta(:,1,:))'); shading flat; colorbar
  title('theta [K]');
  subplot(325); pcolor(squeeze(rh(:,1,:))'); shading flat; colorbar
  title('rh');

  % Use a reference isotopic composition related to vapor in
  %  equilibrium with ocean ater at Tref = 20 deg C
  %   qv =0.8*qsat(p = 1e5 Pa, T = 20 deg C)
  %   Rv = Rocn/alpha(T = 20 deg C)
  Tref = 20 + 273.15;
  RHref = 0.8;
  Rocn = 1; % isotope ratio of ocean water relative to SMOW.
  pref = 1e5;
  theta_ref = Tref*(1e5./pref).^(Rd/Cp);

  q0 = RHref*qsw(1e5,Tref);
  R0_HDO = 1/alpha_equil_liq(Tref,19);
  R0_H2O18 = 1/alpha_equil_liq(Tref,20);
  dD0 = 1000*(R0_HDO-1);
  dO180 = 1000*(R0_H2O18-1);

  [p_Rayleigh,theta_Rayleigh,qv_Rayleigh,delHDO_Rayleigh,delO18_Rayleigh] = ...
      rayleigh_curve(theta_ref,pref,q0,dD0,dO180);

% Start the Rayleigh curve where qv starts decreasing and end it at
% the minimum value of qv in the model domain.
  i_start = min(find(diff(qv_Rayleigh)<0));
  i_end = min(find(qv_Rayleigh<min(qv(:))));
  ind2 = [i_start:i_end];
  qv_Rayleigh = qv_Rayleigh(ind2);
  delHDO_Rayleigh = delHDO_Rayleigh(ind2);
  delO18_Rayleigh = delO18_Rayleigh(ind2);


%bloss: Option to limit depletion in the upper troposphere and stratosphere.
% $$$   % limit depletion in deltaD to -650 per mil.  Find equivalent
% $$$   % value of delO18.
% $$$   dD_min = -650;
% $$$   ind3 = min(find(delHDO_Rayleigh<dD_min));
% $$$   dO18_min = interp1(delHDO_Rayleigh(ind3-1:ind3), ...
% $$$                      delO18_Rayleigh(ind3-1:ind3),dD_min);
% $$$   delHDO_Rayleigh = max(dD_min,delHDO_Rayleigh);
% $$$   delO18_Rayleigh = max(dO18_min,delO18_Rayleigh);
  
  figure(2); clf
  subplot(121); semilogy(delHDO_Rayleigh,qv_Rayleigh)
  subplot(122); plot(delHDO_Rayleigh,delO18_Rayleigh,delHDO_Rayleigh,delHDO_Rayleigh/8)

  dD_vapor = interp1(qv_Rayleigh,delHDO_Rayleigh,qv(:));
  dD_vapor(isnan(dD_vapor)) = -1000;
  hdo_qvapor = (1 + dD_vapor(:)/1000).*qv(:);
  dO18_vapor = interp1(qv_Rayleigh,delO18_Rayleigh,qv(:));
  dO18_vapor(isnan(dO18_vapor)) = -1000;
  o18_qvapor = (1 + dO18_vapor(:)/1000).*qv(:);
  
  whvap = {'dD_','hdo_q','dO18_','o18_q'};
  for m = 1:length(whvap)
    eval(sprintf('%svapor = reshape(%svapor,size(qv));', ...
                whvap{m},whvap{m}));
  end
  ncwrite(in(1).nc,'HDO_QVAPOR',hdo_qvapor);
  ncwrite(in(1).nc,'O18_QVAPOR',o18_qvapor);

  whliq = {'CLOUD','RAIN'};
  for m = 1:length(whliq)
    q = ncread(in(1).nc,['Q' whliq{m}]);
    hdo_q = zeros(size(q));
    o18_q = zeros(size(q));

    ind = find(q>0);
    if ~isempty(ind)
      disp(sprintf('There is some %s in the initial profile',whliq{m}))
      hdo_q(ind) = q(ind).*(1+dD_vapor(ind)/1000).*alpha_equil_liq(T(ind),19);
      o18_q(ind) = q(ind).*(1+dO18_vapor(ind)/1000).*alpha_equil_liq(T(ind),20);
    end
    ncwrite(in(1).nc,['HDO_Q' whliq{m}],hdo_q);
    ncwrite(in(1).nc,['O18_Q' whliq{m}],o18_q);
  end

  whice = {'ICE','SNOW','GRAUP'};
  for m = 1:length(whice)
    q = ncread(in(1).nc,['Q' whice{m}]);
    hdo_q = zeros(size(q));
    o18_q = zeros(size(q));

    ind = find(q>0);
    if ~isempty(ind)
      disp(sprintf('There is some %s in the initial profile',whice{m}))
      hdo_q(ind) = q(ind).*(1+dD_vapor(ind)/1000).*alpha_equil_ice(T(ind),19);
      o18_q(ind) = q(ind).*(1+dO18_vapor(ind)/1000).*alpha_equil_ice(T(ind),20);
    end
    ncwrite(in(1).nc,['HDO_Q' whice{m}],hdo_q);
    ncwrite(in(1).nc,['O18_Q' whice{m}],o18_q);
  end


  figure(3); clf
  subplot(321); pcolor(squeeze(qv(:,1,:))'); shading flat; colorbar
  title('qv [kg kg^{-1}]');
  subplot(322); pcolor(squeeze(T(:,1,:))'); shading flat; colorbar
  title('T [K]');
  subplot(323); pcolor(squeeze(hdo_qvapor(:,1,:))'); shading flat; colorbar
  title('hdo qv [kg kg^{-1}]');
  subplot(324); pcolor(squeeze(dD_vapor(:,1,:))'); shading flat; colorbar
  title('dD vap [per mil]');
  subplot(325); pcolor(squeeze(o18_qvapor(:,1,:))'); shading flat; colorbar
  title('o18 qv [kg kg^{-1}]');
  subplot(326); pcolor(squeeze(dO18_vapor(:,1,:))'); shading flat; colorbar
  title('dO18 vap [per mil]');


otherwise
 error(['Bad option in DuplicateWaterFields.m, choose either 1 or ' ...
        '2.  Stopping']); 


end

error('stop here');

wh = {'VAPOR','CLOUD','ICE','RAIN','SNOW','GRAUP'};
whiso = {'HDO','O18'};
for m = 1:length(wh)
  varname = sprintf('Q%s',wh{m})
  qschema = ncinfo(in(1).nc,varname);
  for nn = 1:length(whiso)
    isoname = sprintf('%s_Q%s',whiso{nn},wh{m})
    qschema.Name = isoname
    for k = 1:length(qschema.Attributes)
      attname = qschema.Attributes(k).Name;
      if strcmp(attname,'description') | ...
        strcmp(attname,'long_name') 
        qschema.Attributes(k).Value = [whiso{nn} ...
                            ' ' qschema.Attributes(k).Value];
      end
    end
    ncwriteschema(in(1).nc,qschema);
  end
end
