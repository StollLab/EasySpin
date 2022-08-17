function ok = test(opt)

% ISC triplet: compare results for population vector provided in xyz or
% energy-ordered basis

D = [900 -600]; % MHz
E = [-280 50]; % MHz

for i = 1:numel(D)

  Triplet.S = 1;
  Triplet.D = [D(i) E(i)];  % MHz
  Triplet.lwpp = 1;  % mT

  Exp.mwFreq = 9.5;  % GHz
  Exp.Range = [290 390];  % mT
  Exp.Harmonic = 0;

  % Provide population vector as [px py pz]
  xyzpops = rand(1,3);
  xyzpops = xyzpops/sum(xyzpops);
  Triplet.initState = {xyzpops,'xyz'};
  [~,spc_xyz] = pepper(Triplet,Exp);

  % Corresponding energy-ordered population vector
  px = xyzpops(1); py = xyzpops(2); pz = xyzpops(3);
  if D(i)>0 && E(i)<0
    p1 = pz; p2 = py; p3 = px;
  elseif D(i)>0 && E(i)>0
    p1 = pz; p2 = px; p3 = py;
  elseif D(i)<0 && E(i)>0
    p1 = px; p2 = py; p3 = pz;
  elseif D(i)<0 && E(i)<0
    p1 = py; p2 = px; p3 = pz;
  end
  Triplet.initState = {[p1 p2 p3],'zerofield'};
  [B,spc_Eorder] = pepper(Triplet,Exp);

  if opt.Display
    hold on
    plot(B,spc_xyz,B,spc_Eorder)
  end

  ok(i) = areequal(spc_xyz,spc_Eorder,1e-16,'rel');

end

end
