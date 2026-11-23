function ok = test() %#ok<FNDEF>

% Field-reversed orientations on partial-sphere powder grids
% - full engine: with a nuclear quadrupole interaction or nuclear-nuclear
%   couplings, the signals at the field directions n and -n differ; the
%   partial grid (Ci) with the added antipodal orientations gives the same
%   signal as the full sphere (C1): 14N 3p ESEEM, 2H Mims ENDOR, 19F pair
%   with nn
% - full engine, initial states: a time-odd electron-only state cancels
%   between n and -n, so its time-even part replaces the orientations -n
%   (S = 1 echo from a state with a polarization along a molecular axis;
%   zero signal for S = 1/2); for two electrons with an antisymmetric
%   correlation, the orientations -n are added
% - full engine: no partners where n and -n give the same signal (1H,
%   eigenstate populations)
% - both engines on the default grid: 14N 3p ESEEM

GridOpt = struct('Method','full','GridSize',7,'GridFrame',[0 0 0]);

p90 = struct('Flip',pi/2,'tp',0.016);  % µs
i90 = struct('Flip',pi/2);  % ideal pulses
i180 = struct('Flip',pi);

% 14N with Q, 3p ESEEM with finite pulses
Sys = struct('Nucs','14N','A',[1.6 1.9 2.6],'AFrame',[0.3 0.5 0.2],'Q',[0.9 0.6]);
Exp.Field = 324.9;  % mT
Exp.mwFreq = mt2mhz(Exp.Field,gfree)/1e3;  % GHz
Exp.Sequence = {p90 0.12 p90 0.1 p90 0.12};
Exp.Dim1 = {'d2' 0.05};  % µs
Exp.nPoints = 4;
Exp.Pathways = '+0-';
ok(1) = partialEqualsFull(Sys,Exp,GridOpt,1e-10);

% 2H with Q, Mims ENDOR with an RF pi pulse, at W band (NMR frequencies
% well above the nutation frequency: RF pulses in the RF frame)
Sys = struct('Nucs','2H','A',[-0.4 -0.4 0.9],'AFrame',[0 0.6 0.3],'Q',[0.2 0.1]);
Exp = struct('Field',3390,'mwFreq',95);
Exp.Sequence = {p90 0.3 p90 0 struct('Channel','rf','tp',4,'Flip',pi) 0.1 p90 0.3};
Exp.rfFreq = 21.9;  % MHz
Exp.Dim1 = {'rfFreq',0.2};
Exp.nPoints = 3;
Exp.Pathways = '+00-';
ok(2) = partialEqualsFull(Sys,Exp,GridOpt,1e-10);

% two 19F with nn, 2p ESEEM with ideal pulses
Sys = struct('Nucs','19F,19F','A',[1 2 5; 0.5 -1 3],'AFrame',[0 0 0; 0.4 0.9 0.2],...
  'nn',[0.3 0.4 -0.7],'nnFrame',[0.2 0.6 0]);
Exp = struct('Field',350,'mwFreq',9.81);
Exp.Sequence = {i90 0.12 i180 0.12};
Exp.Dim1 = {'d1,d2' 0.05};  % µs
Exp.nPoints = 4;
Exp.Pathways = '+-';
ok(3) = partialEqualsFull(Sys,Exp,GridOpt,1e-10);

% S = 1, state with a polarization along a molecular axis (time-odd) and
% an alignment (time-even), 2p echo with finite pulses: without nuclei, the
% time-odd part cancels between n and -n, and the time-even part of the
% state replaces the orientations -n
Sys = struct('S',1,'D',[300 50]);
Sys.initState = {diag([0.6 0.3 0.1]) + 0.1*[0 -1i 0; 1i 0 -1i; 0 1i 0],'uncoupled'};
a90 = struct('tp',0.02,'Amplitude',6);  % µs, MHz
a180 = struct('tp',0.04,'Amplitude',6);
Exp = struct('Field',345,'mwFreq',9.5);
Exp.Sequence = {a90 0.15 a180 0.15};
Exp.Dim1 = {'d1,d2' 0.05};  % µs
Exp.nPoints = 4;
ok(4) = partialEqualsFull(Sys,Exp,GridOpt,1e-10,"none");

% no partners: 1H without Q, with -Sz and with eigenstate populations,
% 3p ESEEM with ideal pulses
Sys = struct('g',[2.0 2.02 2.04],'Nucs','1H','A',[2 6],'AFrame',[0.5 0.3 0.7]);
Exp = struct('Field',345,'mwFreq',9.81);
Exp.Sequence = {i90 0.12 i90 0.1 i90 0.12};
Exp.Dim1 = {'d2' 0.05};  % µs
Exp.nPoints = 4;
Exp.Pathways = '+0-';
ok(5) = partialEqualsFull(Sys,Exp,GridOpt,1e-6,"none");
Sys.initState = {[0.8 0.2],'eigen'};
ok(6) = partialEqualsFull(Sys,Exp,GridOpt,1e-6,"none");
% polarization along a molecular axis, the only electron state of S = 1/2
% besides the identity: time-odd, so the powder signal vanishes
Sys.initState = {-[1 0; 0 -1]/2,'uncoupled'};
[~,y,info] = saffron(Sys,Exp,setfield(GridOpt,'GridSymmetry','Ci')); %#ok<SFLD>
ok(8) = info.reversedField=="none" && all(abs(y)<1e-12);

% two electrons, time-even state with an antisymmetric correlation
% (S1 x S2)_y: its flip-flop elements in the lab frame are imaginary, so
% the orientations -n are needed (without them, the partial grid gives a
% signal 500 times too large). The powder signal is small (4e-4 of that of
% -Sz), so its relative rounding errors are larger.
Sys = struct('S',[1/2 1/2],'g',[2.0 2.01 2.02; 2.003 2.006 2.03],...
  'gFrame',[0 0 0; 0.3 0.6 0.2],'J',10,'dip',[-10 -10 20]);
Sys.initState = {eye(4)/4 + 0.4*(sop([1/2 1/2],'zx')-sop([1/2 1/2],'xz')),'uncoupled'};
Exp = struct('Field',339,'mwFreq',9.5);
b90 = struct('tp',0.01,'Amplitude',25);  % µs, MHz
b180 = struct('tp',0.02,'Amplitude',25);
Exp.Sequence = {b90 0.15 b180 0.15};
Exp.Dim1 = {'d1,d2' 0.05};  % µs
Exp.nPoints = 4;
ok(9) = partialEqualsFull(Sys,Exp,GridOpt,1e-7);

% both engines on the default grid (D2h for this system): with
% field-reversed partners, they agree within the high-field approximation
% and the time step of the full engine, as on a single orientation
Sys = struct('Nucs','14N','A',[1.6 1.9 2.6],'AFrame',[0.3 0.5 0.2],'Q',[0.9 0.6]);
Exp = struct('Field',324.9);  % mT
Exp.mwFreq = mt2mhz(Exp.Field,gfree)/1e3;  % GHz
Exp.Sequence = {p90 0.12 p90 0.1 p90 0.12};
Exp.Dim1 = {'d2' 0.05};  % µs
Exp.nPoints = 6;
Exp.Pathways = '+0-';
O = struct('GridSize',7);
O.Method = 'manifold';
[~,y1,info1] = saffron(Sys,Exp,O);
O.Method = 'full';
[~,y2,info2] = saffron(Sys,Exp,O);
ok(7) = info1.reversedField=="compute" && info2.reversedField=="compute" && ...
  areequal(y1,y2,1e-4,'rel');

end

% partial-sphere grid (Ci) against the full sphere (C1), and the handling
% of the field-reversed partners
function ok = partialEqualsFull(Sys,Exp,Opt,tol,reversedField)
if nargin<5, reversedField = "compute"; end
Opt.GridSymmetry = 'Ci';
[~,yPartial,info] = saffron(Sys,Exp,Opt);
Opt.GridSymmetry = 'C1';
[~,yFull] = saffron(Sys,Exp,Opt);
ok = info.reversedField==reversedField && areequal(yPartial,yFull,tol,'rel');
end
