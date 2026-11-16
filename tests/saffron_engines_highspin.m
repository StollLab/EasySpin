function ok = test() %#ok<FNDEF>

% Electron spins S > 1/2
% - microwave pulses need Amplitude (nutation frequency nu1 in MHz, pulse
%   Hamiltonian nu1*Sx); Flip, Qcrit, ideal pulses, flip-angle increments,
%   predefined experiments and Opt.IdealPulses give errors in both engines
% - with selective rectangular pulses, transition a<->b nutates with
%   nu1*|<a|S+|b>|: the subspace and full-space engines agree on the
%   echoes of the allowed transitions of S = 1 and S = 3/2
% - signal scale of the subspace engine with a core system (S = 1/2 and a
%   nucleus outside the nuclear subspace): primary echo 1/2

% Errors
SysE.S = 1;
SysE.D = 300;  % MHz
SysE.Nucs = '1H';
SysE.A = [0.5 0.5 1];  % MHz
Exp.Field = 340;  % mT
Exp.mwFreq = 9.5;  % GHz
Exp.MolFrame = [0 0.3 0];
Exp.CrystalSymmetry = 1;
rect.tp = 0.06;  % µs
rect.Amplitude = 2;  % MHz
methods = {struct, struct('Method','subspace'), struct('Method','fullspace')};
for m = 1:numel(methods)
  O = methods{m};
  E = Exp; p = rect; p.Flip = pi/2; E.Sequence = {p 0.2};
  ok(m,1) = throwsWith(@() saffron(SysE,E,O),'Amplitude'); %#ok<AGROW>
  E = Exp; p = rmfield(rect,'Amplitude'); p.Flip = pi/2; E.Sequence = {p 0.2};
  ok(m,2) = throwsWith(@() saffron(SysE,E,O),'Amplitude'); %#ok<AGROW>
  E = Exp; p = rmfield(rect,'Amplitude'); p.Type = 'quartersin/linear'; p.trise = 0.01;
  p.Frequency = [-50 50]; p.Qcrit = 5; E.Sequence = {p 0.2};
  ok(m,3) = throwsWith(@() saffron(SysE,E,O),'Qcrit'); %#ok<AGROW>
  E = Exp; p = rmfield(rect,'tp'); E.Sequence = {p 0.2};
  ok(m,4) = throwsWith(@() saffron(SysE,E,O),'ideal pulses'); %#ok<AGROW>
  E = Exp; E.Sequence = {rect 0.2 rect 0.2}; E.Dim1 = {'p2.Flip' 0.1}; E.nPoints = 3;
  ok(m,5) = throwsWith(@() saffron(SysE,E,O),'p2.Amplitude'); %#ok<AGROW>
  E = Exp; E.Sequence = '2pESEEM'; E.dt = 0.01; E.nPoints = 4;
  ok(m,6) = throwsWith(@() saffron(SysE,E,O),'Predefined experiments'); %#ok<AGROW>
  E = Exp; E.Sequence = {rect 0.2 rect 0.2};
  ok(m,7) = throwsWith(@() saffron(SysE,E,setfield(O,'IdealPulses',true)),'IdealPulses'); %#ok<AGROW>
end
ok = ok(:).';

% Engine agreement on the allowed transitions, primary echo with selective
% rectangular pulses (two-level systems for the subspace engine)
ph = (0:3).'*pi/2;
sel = [ph 0*ph -ph];  % echo pathway by phase cycle: first pulse in steps of 90°
a90.tp = 0.06; a90.Amplitude = 2;  % µs, MHz
a180.tp = 0.12; a180.Amplitude = 2;
Exp.Sequence = {a90 0.15 a180 0.15};
Exp.Dim1 = {'d1,d2' 0.02};
Exp.nPoints = 3;
Exp.DetPhase = -pi/2;
for S = [1 3/2]
  Sys = SysE;
  Sys.S = S;
  ExpCW = struct('mwFreq',Exp.mwFreq,'Range',[250 450],'MolFrame',Exp.MolFrame,'CrystalSymmetry',1);
  B = resfields(Sys,ExpCW);  % allowed transitions (Delta mS = 1)
  for iB = 1:numel(B)
    E = Exp;
    E.Field = B(iB);
    Es = E; Es.Pathways = 'echo';
    [~,ySub] = saffron(Sys,Es,struct('Method','subspace'));
    Ef = E; Ef.PhaseCycle = sel; Ef.DetWindow = 0;
    % simulation frame close to the microwave frequency: long time steps
    [~,yFull] = saffron(Sys,Ef,struct('Method','fullspace','SimFreq',9.4));
    % high-field approximation of the subspace engine, ~(D/nu_e)^2
    ok(end+1) = areequal(ySub,yFull,1e-2,'rel'); %#ok<AGROW>
  end
end

% The transition moment matters much more than the tolerance: the same
% echo with nu1 scaled by 1/sqrt(2) (S = 1, without the moment)
Sys = SysE;
B = resfields(Sys,ExpCW);
E = Exp; E.Field = B(1); E.Pathways = 'echo';
[~,y1] = saffron(Sys,E,struct('Method','subspace'));
E.Sequence{1}.Amplitude = 2/sqrt(2);
E.Sequence{3}.Amplitude = 2/sqrt(2);
[~,y2] = saffron(Sys,E,struct('Method','subspace'));
ok(end+1) = max(abs(y1-y2))/max(abs(y1)) > 0.1;

% S = 1/2 with a nucleus in the core system (Opt.Nuclei): primary echo 1/2,
% up to the hyperfine mixing in the core system, (A/nu_e)^2 ~ 1e-8
Sys12.g = gfree;
Sys12.Nucs = '1H,1H';
Sys12.A = [0 0 0; 2 2 2];  % MHz
p90.Flip = pi/2;
p180.Flip = pi;
E = Exp;
E.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
E.Sequence = {p90 0.15 p180 0.15};
E.Pathways = 'echo';
[~,y] = saffron(Sys12,E,struct('Method','subspace','Nuclei',1));
ok(end+1) = areequal(y,0.5*ones(size(y)),1e-6,'abs');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch ME
  ok = contains(ME.message,str);
end
end
