function [err,data] = test(opt,olddata)

% S=3/2 regression test

Sys = struct('S',3/2,'g',[2 2 2],'D',[-1 -1 2]*3e3);
Exp = struct('mwFreq',9.5,'Range',[0 2000]);
[phi,theta] = sphgrid('Dinfh',101);
z = zeros(size(phi));
Exp.CrystalOrientation = [phi(:) theta(:) z(:)];
[p,i] = resfields(Sys,Exp);

if (opt.Display)
  title('Test 3: S=3/2 system');
  plot(theta,p,'.');
end

data.p = p;
data.i = i;

if isempty(olddata)
  err = [];
else
  idxp = ~isnan(p);
  idxi = ~isnan(i);  
  ok = areequal(p(idxp),data.p(idxp),1e-4,'abs') & areequal(i(idxi),data.i(idxi),1e-4,'abs');
  err = ~ok;
end
