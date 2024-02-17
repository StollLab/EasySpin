function ok = diptensor_en()

% Test whether diptensor gives correct results for hyperfine coupling tensors.

Isotopes = {'14N','1H','15N','13C','2H'};
rvec = [0;0;1]*0.2; % nm
r = norm(rvec)*1e-9; % nm -> m

for k = 1:numel(Isotopes)
  
  T = diptensor(gfree,Isotopes{k},rvec);
  Tperp = T(1,1);
  
  gn = nucgval(Isotopes{k});
  Tperp0 = -(mu0/4/pi)*r^-3*bmagn*gfree*nmagn*gn;  % J
  Tperp0 = Tperp0/planck/1e6;  % J -> MHz
  
  ok(k) = areequal(Tperp,Tperp0,1e-10,'rel');
    
end
