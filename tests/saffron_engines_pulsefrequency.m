function ok = test() %#ok<FNDEF>

% Finite microwave pulse with a frequency offset f from Exp.mwFreq
% (Pulse.Frequency). The subspace engine computes its propagator in the
% frame of the pulse and transforms it back to the frame rotating at
% Exp.mwFreq.
% - both engines agree, at f = 0 and f ~= 0; without the back-transform,
%   the phase would be off by 2*pi*f*tp
% - a pulse with Flip = 0 is free evolution, independent of f
% - powder with a nuclear quadrupole interaction: the field-reversed
%   partners of a partial-sphere grid give the signal of the full sphere

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [3 8];  % MHz

Exp.Field = 1200;  % mT
Exp.mwFreq = mt2mhz(Exp.Field,Sys.g)/1e3;  % GHz, on resonance
Exp.SampleFrame = [0 -0.7 0];
p90 = struct('Type','rectangular','Flip',pi/2,'tp',0.016);  % µs
p180 = struct('Type','rectangular','Flip',pi,'tp',0.032);  % µs

tol = 1e-4;
fList = [0 5 -8];  % MHz
for iF = 1:numel(fList)
  p90.Frequency = fList(iF);
  Exp.Sequence = {p90 0.3 p180 0.3};
  [~,ySub] = saffron(Sys,Exp,struct('Method','subspace'));
  [~,yFull] = saffron(Sys,Exp,struct('Method','fullspace'));
  phaseError = abs(exp(2i*pi*fList(iF)*p90.tp)-1);  % without back-transform
  ok(iF) = areequal(ySub,yFull,tol,'rel') && (fList(iF)==0 || phaseError>100*tol); %#ok<AGROW>
end

% Flip = 0: free evolution, independent of the frequency offset
p0 = struct('Type','rectangular','Flip',0,'tp',0.032,'Frequency',0);
p90.Frequency = 0;
Exp.Sequence = {p90 0.2 p0 0.1 p180 0.3};
[~,y0] = saffron(Sys,Exp,struct('Method','subspace'));
Exp.Sequence{3}.Frequency = 5;  % MHz
[~,y5] = saffron(Sys,Exp,struct('Method','subspace'));
ok(end+1) = areequal(y5,y0,1e-10,'rel');

% Powder, 14N with Q: partial-sphere grid with field-reversed partners
% (computed explicitly for finite pulses) against the full sphere
SysN.Nucs = '14N';
SysN.A = [1.6 1.9 2.6];  % MHz
SysN.AFrame = [0.3 0.5 0.2];
SysN.Q = [0.9 0.6];  % MHz
ExpN.Field = 324.9;  % mT
ExpN.mwFreq = mt2mhz(ExpN.Field,gfree)/1e3;  % GHz
f90 = struct('Flip',pi/2,'tp',0.016,'Frequency',4);  % µs, MHz
f180 = struct('Flip',pi,'tp',0.032);  % µs
ExpN.Sequence = {f90 0.1 f180 0.1};
ExpN.Dim1 = {'d1,d2' 0.05};  % µs
ExpN.nPoints = 16;
Opt = struct('Method','subspace','GridSize',7,'GridFrame',[0 0 0]);
Opt.GridSymmetry = 'Ci';
[~,yPartial,info] = saffron(SysN,ExpN,Opt);
Opt.GridSymmetry = 'C1';
[~,yFull] = saffron(SysN,ExpN,Opt);
ok(end+1) = info.reversedField=="compute" && areequal(yPartial,yFull,1e-10,'rel');
