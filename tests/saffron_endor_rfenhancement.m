function ok = test() %#ok<FNDEF>

% Hyperfine enhancement of RF pulses (Opt.rfEnhancement): the RF field
% also tilts the electron spin, which changes the hyperfine field at the
% nucleus. By default, both engines include this to first order in 1/nu_e;
% with Opt.rfEnhancement = false, the RF field acts only through the
% nuclear Zeeman interaction.
% - A flip angle (Flip) sets the nominal field B2 for the first nucleus,
%   so it gives the same signal as that B2, with and without enhancement.
% - The enhancement changes the ENDOR lines substantially for 1H with A
%   comparable to nuI.
% - Full engine: the first-order operator in the shifted frame agrees with
%   the lab frame, where the RF field acts on the full magnetic moment
%   (exact).
% - Manifold engine with nuclei in the core system (here 14N, Opt.Nuclei):
%   the enhancement from the tilt of the core eigenstates agrees with the
%   two-manifold model.
% - Natural-abundance nuclei (Sys.Nucs = 'H') define the nominal B2 by
%   their most abundant isotope.

p90.Flip = pi/2;
Sys.Nucs = '1H';
Sys.A = [4 4 10];  % MHz
rf.Channel = 'rf';
rf.tp = 2;  % µs
rf.Flip = pi;
Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0.6 0];
Exp.CrystalSymmetry = 1;
Exp.Sequence = {p90 0.2 p90 0 rf 0.1 p90 0.2};
Exp.Pathways = '+00-';
Exp.DetPhase = -pi/2;
% RF frequencies at the centers of the two ENDOR lines, nu(mS) =
% |mS*A*z - nuI*z| with the field direction z in the molecular frame
nuI = larmorfrq('1H',Exp.Field);
z = [sin(0.6); 0; cos(0.6)];
nuLine = @(mS) norm(mS*diag(Sys.A)*z - nuI*z);  % MHz
Exp.rfFreq = nuLine(+1/2);
Exp.nPoints = 2;
Exp.Dim1 = {'rfFreq',nuLine(-1/2)-nuLine(+1/2)};

% Flip and the equivalent B2
rfB = rmfield(rf,'Flip');
rfB.B2 = rf.Flip/(2*pi*rf.tp*larmorfrq('1H',1));  % mT
ExpB = Exp;
ExpB.Sequence{5} = rfB;
methods = {'manifold','full'};
ok = false(1,0);
y = cell(2,2);
for m = 1:2
  for enh = [true false]
    Opt = struct('Method',methods{m},'rfEnhancement',enh);
    [~,yFlip] = saffron(Sys,Exp,Opt);
    [~,yB2] = saffron(Sys,ExpB,Opt);
    ok(end+1) = areequal(yFlip,yB2,1e-8,'rel'); %#ok<AGROW>
    y{m,2-enh} = yFlip;
  end
  % the enhancement changes the signal by much more than the tolerances below
  ok(end+1) = max(abs(y{m,1}-y{m,2}))>0.01*max(abs(y{m,1})); %#ok<AGROW>
end

% Full engine: shifted frame against the lab frame, change of the signal by
% the enhancement. At the line centers, the rotating-wave error of the
% shifted frame (Bloch-Siegert shift of the lines) changes the signal only
% to second order, and the frames agree to about 1% of the change.
ExpL = Exp;
ExpL.rfPhaseLocked = true;  % one RF phase, keeps the lab-frame runs short
dy = cell(1,2);
for lab = [false true]
  yE = cell(1,2);
  for enh = [true false]
    Opt = struct('Method','full','rfEnhancement',enh);
    if lab, Opt.SimFreq = 0; end
    [~,yE{2-enh}] = saffron(Sys,ExpL,Opt);
  end
  dy{1+lab} = yE{1} - yE{2};
end
ok(end+1) = areequal(dy{1},dy{2},0.02*max(abs(dy{2})),'abs');

% Manifold engine, 14N in the core system: the 1H ENDOR effect relative to
% the echo agrees with the spin system without 14N
SysN = Sys;
SysN.Nucs = '1H,14N';
SysN.A = [4 4 10; 20 20 60];  % MHz
OptN = struct('Method','manifold','Nuclei',1);
Exp0 = Exp;
Exp0.Sequence{5}.Flip = 0;
for enh = [true false]
  OptN.rfEnhancement = enh;
  [~,yN] = saffron(SysN,Exp,OptN);
  [~,yN0] = saffron(SysN,Exp0,OptN);
  Opt = struct('Method','manifold','rfEnhancement',enh);
  [~,yH] = saffron(Sys,Exp,Opt);
  [~,yH0] = saffron(Sys,Exp0,Opt);
  ok(end+1) = areequal(yN/yN0(1)-1,yH/yH0(1)-1,2e-4,'abs'); %#ok<AGROW>
end

% Natural-abundance nuclei
SysH = Sys;
SysH.Nucs = 'H';
for m = 1:2
  [~,yH] = saffron(SysH,Exp,struct('Method',methods{m}));
  ok(end+1) = areequal(yH,y{m,1},1e-3,'rel'); %#ok<AGROW>
end

% Opt.rfEnhancement must be true or false
try
  saffron(Sys,Exp,struct('rfEnhancement',2));
  ok(end+1) = false;
catch err
  ok(end+1) = contains(err.message,'Opt.rfEnhancement must be true or false');
end
