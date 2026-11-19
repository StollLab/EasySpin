function ok = test() %#ok<FNDEF>

% Electron-only initial states (Sys.initState) in the subspace engine
% - powders: the subspace and full-space engines agree, for S = 1/2
%   (eigenstate populations) and for S = 1 (zero-field triplet populations,
%   selective pulses)
% - a time-odd state (polarization along a molecular axis) averages to
%   zero over the sphere; a partial-sphere grid (Ci) gives the same signals
%   as the full sphere (C1), also with field-reversed partners (Q)
% - states with a nuclear part, crystals, and relaxation need the
%   full-space engine

GridOpt = struct('GridSymmetry','C1','GridFrame',[0 0 0],'GridSize',7);

% S = 1/2: populations of the eigenstates, ideal pulses
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz
Exp.Field = 350;  % mT
Exp.mwFreq = mt2mhz(Exp.Field,Sys.g)/1e3;  % GHz
p90.Flip = pi/2;
p180.Flip = pi;
Exp.Sequence = {p90 0.2 p180 0.2};
Exp.DetPhase = -pi/2;

Sys.initState = {[0.8 0.2],'eigen'};  % lower and upper level
[~,ySub,info] = saffron(Sys,Exp,setfield(GridOpt,'Method','subspace')); %#ok<SFLD>
[~,yFull] = saffron(Sys,Exp,setfield(GridOpt,'Method','fullspace')); %#ok<SFLD>
[~,yDefault] = saffron(rmfield(Sys,'initState'),Exp,GridOpt);
ok(1) = strcmp(info.Method,'subspace') && areequal(ySub,yFull,1e-5,'rel');
% polarization 0.6 instead of 1, and the trace divided by the number of
% nuclear states
ok(2) = areequal(ySub,yDefault*0.6/2,1e-10,'rel');

% S = 1/2, time-odd state: zero powder signal in both engines and on both
% grids. On the partial-sphere grid without partners, the weights of n and
% -n cancel, so the entries are dropped (one is kept).
Sys.initState = {-[1 0; 0 -1]/2,'uncoupled'};  % -Sz along the molecular z axis
[~,ySub] = saffron(Sys,Exp,setfield(GridOpt,'Method','subspace')); %#ok<SFLD>
[~,yFull] = saffron(Sys,Exp,setfield(GridOpt,'Method','fullspace')); %#ok<SFLD>
[~,yCi,info] = saffron(Sys,Exp,setfield(GridOpt,'GridSymmetry','Ci')); %#ok<SFLD>
ok(3) = abs(ySub)<1e-10*abs(yDefault) && abs(yFull)<1e-6*abs(yDefault) && abs(yCi)<1e-10*abs(yDefault) && ...
  info.nEntries==1;

% S = 1, zero-field populations of the triplet, selective pulses at an
% allowed transition (two-level systems in the subspace engine). The
% engines differ by the high-field approximation: ~4e-3 for -Sz, ~7e-3
% here, since the subspace engine uses the populations of the eigenstates
% (mS mixed to first order in D/nu_e) and the full-space engine removes
% the coherences between different mS (secular eigenstates)
SysT.S = 1;
SysT.D = [300 50];  % MHz
SysT.Nucs = '1H';
SysT.A = [0.5 1];  % MHz
SysT.initState = {[0.5 0.3 0.2],'xyz'};
ExpT.mwFreq = 9.5;  % GHz
ExpT.Field = 345;  % mT
a90 = struct('tp',0.06,'Amplitude',2);  % µs, MHz
a180 = struct('tp',0.12,'Amplitude',2);
ExpT.Sequence = {a90 0.15 a180 0.15};
ExpT.DetPhase = -pi/2;
Es = ExpT; Es.Pathways = 'echo';
[~,ySub] = saffron(SysT,Es,setfield(GridOpt,'Method','subspace')); %#ok<SFLD>
ph = (0:3).'*pi/2;
Ef = ExpT; Ef.PhaseCycle = [ph 0*ph -ph];  % echo pathway
[~,yFull] = saffron(SysT,Ef,struct('Method','fullspace','SimFreq',9.4,...
  'GridSymmetry','C1','GridFrame',[0 0 0],'GridSize',7));
ok(4) = areequal(ySub,yFull,1e-2,'rel') && abs(ySub)>0;

% Partial-sphere grid with field-reversed partners (14N with Q): the same
% signal as the full sphere, for a time-even state (eigenstates: equal
% weights, so the shortcut of the ideal 2pESEEM applies) and a time-odd
% state (partners computed with their own weights), also for S = 1 (core
% system; finite pulses, so the partners are always computed)
SysN.g = 2;
SysN.Nucs = '14N';
SysN.A = [1.6 1.9 2.6];  % MHz
SysN.AFrame = [0.3 0.5 0.2];
SysN.Q = [0.9 0.6];  % MHz
ExpN = Exp;
ExpN.Field = 324.9;  % mT
ExpN.mwFreq = mt2mhz(ExpN.Field,SysN.g)/1e3;
ExpN.Sequence = {p90 0.15 p180 0.15};
ExpN.Dim1 = {'d1,d2' 0.05};  % µs
ExpN.nPoints = 8;
SysT2 = SysT;
SysT2.Nucs = '14N';
SysT2.A = SysN.A; SysT2.AFrame = SysN.AFrame; SysT2.Q = SysN.Q;
ExpT2 = ExpT;
ExpT2.Pathways = 'echo';
ExpT2.Dim1 = {'d1,d2' 0.05};  % µs
ExpT2.nPoints = 8;
cases = {
  SysN, {[0.8 0.2],'eigen'}, ExpN, "conj"
  SysN, {[0.3 0.1-0.2i; 0.1+0.2i -0.3],'uncoupled'}, ExpN, "compute"
  SysT2, {[0.5 0.3 0.2],'xyz'}, ExpT2, "compute"
  SysT2, {diag([0.6 0.3 0.1]) + 0.1*[0 -1i 0; 1i 0 -1i; 0 1i 0],'uncoupled'}, ExpT2, "compute"
  };
for k = 1:size(cases,1)
  S_ = cases{k,1};
  S_.initState = cases{k,2};
  O = GridOpt; O.Method = 'subspace';
  [~,yFull] = saffron(S_,cases{k,3},O);
  O.GridSymmetry = 'Ci';
  [~,yPartial,info] = saffron(S_,cases{k,3},O);
  ok(end+1) = info.reversedField==cases{k,4} && areequal(yPartial,yFull,1e-10,'rel'); %#ok<AGROW>
end

% No excited orientation (excitation window far from resonance): no signal
S_ = SysN;
S_.initState = {[0.8 0.2],'eigen'};
E_ = ExpN;
E_.mwFreq = ExpN.mwFreq + 0.5;  % GHz
E_.Sequence{1}.tp = 0.02;  % µs, excitation window 50 MHz
E_.Sequence{3}.tp = 0.02;
O = GridOpt; O.Method = 'subspace'; O.GridSymmetry = 'Ci'; O.IdealPulses = true;
[~,y,info] = saffron(S_,E_,O);
ok(end+1) = all(y==0) && info.nEntries==0;

% Engine selection
Sys.initState = {eye(4)/4 - kron([1 0; 0 -1]/2,eye(2))/2,'uncoupled'};  % with nuclear part
[~,~,info] = saffron(Sys,Exp);
ok(end+1) = strcmp(info.Method,'fullspace') && contains(info.MethodReason,'nuclear part');
Sys.initState = {[0.8 0.2],'eigen'};
Ec = Exp; Ec.MolFrame = [0 0.3 0]; Ec.CrystalSymmetry = 1;
[~,~,info] = saffron(Sys,Ec);
ok(end+1) = strcmp(info.Method,'fullspace') && contains(info.MethodReason,'crystals');
Sys.T2 = 1;  % µs
try
  saffron(Sys,Exp,struct('Method','subspace'));
  ok(end+1) = false;
catch err
  ok(end+1) = contains(err.message,'initial states with relaxation');
end

end
