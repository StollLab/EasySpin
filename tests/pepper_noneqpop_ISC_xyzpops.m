function ok = test(opt)

% ISC triplet: compare results for population vector provided in xyz or
% energy-ordered basis

rng(6445)

D = 900*[1  1 -1 -1];  % MHz
E = 280*[1 -1  1 -1];  % MHz

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [290 390];  % mT
Exp.Harmonic = 0;

Triplet.S = 1;
Triplet.lwpp = 1;  % mT
Triplet.DFrame = rand(1,3);

xyzpops = rand(1,3);

for i = 1:numel(D)

  Triplet.D = [D(i) E(i)];  % MHz

  % Provide populations for Tx, Ty, Tz
  Triplet.initState = {xyzpops,'xyz'};
  [B,spc_xyz] = pepper(Triplet,Exp);

  % Get energy-ordered population vector
  if D(i)>0
    if E(i)<0
      energyorder = [3 2 1];
    else
      energyorder = [3 1 2];
    end
  else
    if E(i)>0
      energyorder = [1 2 3];
    else
      energyorder = [2 1 3];
    end
  end
  zfpops = xyzpops(energyorder);

  Triplet.initState = {zfpops,'zerofield'};
  [~,spc_zf] = pepper(Triplet,Exp);

  if opt.Display
    hold on
    plot(B,spc_xyz,B,spc_zf,'--')
    title(sprintf('D = %g MHz, E= %g MHz',D(i),E(i)));
    pause
  end

  ok(i) = areequal(spc_xyz,spc_zf,1e-12,'abs');

end

end
