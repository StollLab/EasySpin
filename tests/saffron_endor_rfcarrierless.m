function ok = test() %#ok<FNDEF>

% RF pulses follow the same frequency rule as microwave pulses: the
% frequency is the carrier Exp.rfFreq plus the offset Pulse.Frequency, or,
% without Exp.rfFreq, the absolute Pulse.Frequency. The carrier phase runs
% continuously from the start of the sequence, while the waveform of a
% pulse (with its Frequency) starts at the start of the pulse.

% Manifold engine, Mims ENDOR, single orientation: an RF sweep over the
% carrier, over the offset of the pulse, and over the absolute frequency of
% the pulse (no carrier) give the same spectrum
Sys.Nucs = '1H';
Sys.A = [1 1 3];  % MHz
p90.Flip = pi/2;
rf.Channel = 'rf';
rf.tp = 5;  % µs
rf.Flip = pi;
Exp.Field = 350;  % mT
Exp.MolFrame = [0 0.2 0];
Exp.CrystalSymmetry = 1;
Exp.Sequence = {p90 0.15 p90 0 rf 0 p90 0.15};
Exp.Pathways = '+00-';
Exp.nPoints = 61;
c = 13;    % MHz, carrier
f0 = 0.5;  % MHz, offset of the RF pulse from the carrier
step = 0.05;  % MHz

ExpC = Exp;  % carrier sweep
ExpC.rfFreq = c;
ExpC.Sequence{5}.Frequency = f0;
ExpC.Dim1 = {'rfFreq' step};
[xC,yC,infoC] = saffron(Sys,ExpC);

ExpO = ExpC;  % offset sweep
ExpO.Dim1 = {'p3.Frequency' step};
[xO,yO,infoO] = saffron(Sys,ExpO);

ExpA = rmfield(ExpO,'rfFreq');  % absolute frequency, no carrier
ExpA.Sequence{5}.Frequency = c + f0;
[xA,yA,infoA] = saffron(Sys,ExpA);

ok(1) = strcmp(infoC.Method,'manifold') && strcmp(infoO.Method,'manifold') && ...
  strcmp(infoA.Method,'manifold');
ok(2) = areequal(yO,yC,1e-12,'rel') && areequal(yA,yC,1e-12,'rel') && ...
  any(abs(yC-yC(1))>0.01*abs(yC(1)));
ok(3) = areequal(xC,c+step*(0:Exp.nPoints-1),1e-12,'abs') && ...
  areequal(xO,f0+step*(0:Exp.nPoints-1),1e-12,'abs') && ...
  areequal(xA,c+f0+step*(0:Exp.nPoints-1),1e-12,'abs');

% Full engine: carrier c plus offset and the absolute frequency give
% the same RF field, and the same echo, if the absolute pulse has the phase
% of the carrier at the start of the pulse, 2*pi*c*tStart. (In a single
% orientation, the echo depends on the RF phase if the RF carrier has the
% same phase in every shot.)
ExpF = Exp;
ExpF.mwFreq = 9.8;  % GHz
ExpF.rfPhaseLocked = true;
ExpF.nPoints = 3;
ExpF.Sequence{5} = struct('Channel','rf','tp',2,'B2',0.5/(2*0.04258));  % pi pulse for a bare proton
OptF.Method = 'full';
ExpF1 = ExpF;
ExpF1.rfFreq = c;
ExpF1.Sequence{5}.Frequency = 1.2;
ExpF1.Dim1 = {'p3.Frequency' 0.1};
[~,yF1] = saffron(Sys,ExpF1,OptF);
ExpF2 = rmfield(ExpF1,'rfFreq');
ExpF2.Sequence{5}.Frequency = c + 1.2;
tStart = ExpF.Sequence{2};  % µs, start of the RF pulse (ideal mw pulses)
ExpF2.Sequence{5}.Phase = 2*pi*c*tStart;
[~,yF2] = saffron(Sys,ExpF2,OptF);
ExpF2.Sequence{5}.Phase = 0;
[~,yF2phase0] = saffron(Sys,ExpF2,OptF);
ok(4) = areequal(yF2,yF1,1e-8,'rel') && ~areequal(yF2phase0,yF1,1e-3,'rel');
% with a random RF phase (default), the phase of a single RF pulse does
% not matter
ExpF2.rfPhaseLocked = false;
[~,yF2random] = saffron(Sys,ExpF2,OptF);
ExpF2.Sequence{5}.Phase = 2*pi*c*tStart;
[~,yF2randomPhase] = saffron(Sys,ExpF2,OptF);
okRandomPhase = areequal(yF2randomPhase,yF2random,1e-10,'rel');

% Phase rule, full engine: two RF pi/2 pulses on resonance with one
% NMR transition (A only along z, so the NMR frequencies are nuI -+ A/2),
% acting on nuclear polarization. The other NMR transition is 20 MHz away
% and is not excited. With phases 0 and pi, the second pulse undoes the
% first if the two pulses are coherent at the NMR frequency. The RF
% carrier has the same phase in every shot.
SysP.Nucs = '1H';
SysP.A = [0 0 20];  % MHz
SysP.initState = sop([1/2 1/2],'ez');
ExpP.Field = 350;  % mT
ExpP.mwFreq = 9.8;  % GHz
ExpP.rfPhaseLocked = true;
ExpP.MolFrame = [0 0 0];
ExpP.CrystalSymmetry = 1;
ExpP.DetOperator = {'ez'};
nuRF = larmorfrq('1H',ExpP.Field) + SysP.A(3)/2;  % MHz, resonant NMR transition
rfP.Channel = 'rf';
rfP.tp = 5;  % µs
rfP.Flip = pi/2;
% delay between the pulses: nuRF*dtStart = 1/4 (mod 1), where dtStart is the
% time between the starts of the two pulses
dtStart = (ceil(nuRF*rfP.tp) + 0.25)/nuRF;  % µs
d = dtStart - rfP.tp;
% The final ideal microwave pulse acts on the electron spin only and leaves
% the nuclear polarization unchanged; it gives the sequence detectable
% electron coherence pathways.
seq = @(phi2) {rfP d setfield(rfP,'Phase',phi2) 0 p90 0}; %#ok<SFLD>
echo = @(E) real(saffron(SysP,E,OptF));

E0 = ExpP;  % no RF field: reference
rf0 = rfP; rf0.Flip = 0;
E0.rfFreq = nuRF;
E0.Sequence = {rf0 d rf0 0 p90 0};
y0 = echo(E0);

% carrier at the NMR frequency (offset 0): 0 and pi cancel
Ec = ExpP;
Ec.rfFreq = nuRF;
Ec.Sequence = seq(pi);
yCarrier = echo(Ec);
Ec.Sequence = seq(0);
yCarrier0 = echo(Ec);

% no carrier (absolute frequency): pi does not cancel, pi + 2*pi*f*dtStart does
Ea = ExpP;
Ea.Sequence = seq(pi);
Ea.Sequence{1}.Frequency = nuRF;
Ea.Sequence{3}.Frequency = nuRF;
yAbsPi = echo(Ea);
Ea.Sequence{3}.Phase = pi + 2*pi*nuRF*dtStart;
yAbs = echo(Ea);

% carrier with the pulses at an offset: the same, with the offset frequency
fOff = 0.3;  % MHz
Eo = ExpP;
Eo.rfFreq = nuRF - fOff;
Eo.Sequence = seq(pi);
Eo.Sequence{1}.Frequency = fOff;
Eo.Sequence{3}.Frequency = fOff;
yOffPi = echo(Eo);
Eo.Sequence{3}.Phase = pi + 2*pi*fOff*dtStart;
yOff = echo(Eo);

tolPhase = 1e-3;
ok(5) = areequal(yCarrier,y0,tolPhase,'rel') && abs(yCarrier0-y0)>0.2*abs(y0);
ok(6) = areequal(yAbs,y0,tolPhase,'rel') && abs(yAbsPi-y0)>0.2*abs(y0);
ok(7) = areequal(yOff,y0,tolPhase,'rel') && abs(yOffPi-y0)>0.2*abs(y0);

% RF phase cycle: two steps with phases 0 and pi on the second RF pulse.
% Nuclear polarization (coherence order 0) is detected with the same
% receiver phase in both steps, so the signal is the average of the steps.
Ec.Sequence = seq(0);
Ec.PhaseCycle = [0 0 0 0; 0 pi 0 0];
yCycle = echo(Ec);
ok(8) = areequal(yCycle,(yCarrier0+yCarrier)/2,1e-10,'rel');
ok(9) = okRandomPhase;

% The same sweeps with the values per point (a column), and with the
% carrier sweep split into two lines whose changes add up
ExpO2 = ExpO;
ExpO2.Dim1 = {'p3.Frequency' step*(0:Exp.nPoints-1).'};
[~,yO2] = saffron(Sys,ExpO2);
ExpC2 = ExpC;
ExpC2.Dim1 = {'rfFreq' step/2; 'rfFreq' step/2};
[~,yC2] = saffron(Sys,ExpC2);
ok(end+1) = areequal(yO2,yC,1e-12,'rel') && areequal(yC2,yC,1e-12,'rel');

% Errors
ExpSweepNoCarrier = rmfield(ExpC,'rfFreq');  % carrier sweep without carrier
ExpNoFreq = rmfield(ExpO,'rfFreq');  % neither carrier nor pulse frequency
ExpNoFreq.Sequence{5} = rmfield(ExpNoFreq.Sequence{5},'Frequency');
ExpZero = ExpC;  % zero carrier
ExpZero.rfFreq = 0;
ExpVector = ExpC;  % range in a user-defined sequence
ExpVector.rfFreq = [10 20];
ExpNeg = ExpO;  % 13.5 MHz down to -16.5 MHz
ExpNeg.Dim1 = {'p3.Frequency' -0.5};
ExpChirp = rmfield(ExpO,'rfFreq');  % chirp with band edges moved by -0.5 MHz per step
ExpChirp.Sequence{5}.Frequency = [5 6];
ExpChirp.Dim1 = {'p3.Frequency' [-0.5 -0.5]};
ExpEdge = ExpChirp;  % lower band edge moved by -0.2 MHz per step
ExpEdge.Dim1 = {'p3.Frequency(1)' -0.2};
ExpRow = ExpO;  % a row of values per point is read as band-edge steps
ExpRow.Dim1 = {'p3.Frequency' step*(0:Exp.nPoints-1)};
ExpRange = struct('Field',350,'Sequence','MimsENDOR','tau',0.1,'Range',[10 20]);
ExpRangeESEEM = struct('Field',350,'Sequence','3pESEEM','tau',0.1,'dt',0.01,'Range',[10 20]);
ExpRfESEEM = struct('Field',350,'Sequence','3pESEEM','tau',0.1,'dt',0.01,'rfFreq',[10 20]);
errorCases = {
  ExpSweepNoCarrier, 'Exp.rfFreq (MHz) is needed'
  ExpNoFreq,         'Without Exp.rfFreq'
  ExpZero,           'must be a single positive number'
  ExpVector,         'must be a single positive number'
  ExpNeg,            'must be positive'
  ExpChirp,          'must be positive'
  ExpEdge,           'must be positive'
  ExpRow,            'must be a step'
  ExpRange,          'Exp.Range is no longer supported'
  ExpRangeESEEM,     'Exp.Range is not a parameter'
  ExpRfESEEM,        'Exp.rfFreq cannot be used'
  };
for k = 1:size(errorCases,1)
  try
    saffron(Sys,errorCases{k,1});
    okErr = false;
  catch err
    okErr = contains(err.message,errorCases{k,2});
  end
  ok(end+1) = okErr; %#ok<AGROW>
end
