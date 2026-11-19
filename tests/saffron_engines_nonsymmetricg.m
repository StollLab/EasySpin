function ok = test() %#ok<FNDEF>

% S = 1/2 with a non-symmetric g matrix: the Zeeman term bmagn*B*zL.'*g*S
% quantizes the electron spin along g.'*zL, which sets the EPR frequency
% and the direction of <S> in the hyperfine terms. Check of the EPR
% frequency of the subspace engine against resfreqs_matrix, through the
% excitation window of ideal pulses (Opt.IdealPulses, FWHM 1/tp = 5 MHz):
% on resonance, the echo has the full amplitude, and 2.5 MHz off
% resonance half of it. With g*zL instead, the frequency would be off by
% about 6 MHz.

Sys.g = [2.0 0.1 0.05; -0.08 2.1 0.02; 0.03 -0.04 2.2];
Sys.Nucs = '1H';
Sys.A = [1 1.5 3];  % MHz, below the window width: nucleus in the nuclear subspace
Sys.AFrame = [0.2 0.6 0];

Exp.Field = 350;  % mT
Exp.MolFrame = [0.4 1.0 -0.3];
Exp.CrystalSymmetry = 1;
SysE = rmfield(Sys,{'Nucs','A','AFrame'});
nuEPR = resfreqs_matrix(SysE,Exp);  % MHz

p90 = struct('Flip',pi/2,'tp',0.2);  % µs
p180 = struct('Flip',pi,'tp',0.2);  % µs
Exp.Sequence = {p90 0.3 p180 0.3};
Exp.Pathways = 'echo';  % refocused: independent of the resonance offset
Exp.DetPhase = -pi/2;
Opt.Method = 'subspace';

% ideal pulses without excitation window
E0 = Exp;
E0.Sequence{1}.tp = 0;
E0.Sequence{3}.tp = 0;
E0.mwFreq = nuEPR/1e3;  % GHz
[~,y0] = saffron(Sys,E0,Opt);

Opt.IdealPulses = true;
E = Exp;
E.mwFreq = nuEPR/1e3;  % GHz
[~,yRes] = saffron(Sys,E,Opt);
E.mwFreq = (nuEPR+2.5)/1e3;  % GHz
[~,yOff] = saffron(Sys,E,Opt);

ok(1) = areequal(yRes,y0,1e-6,'rel');
ok(2) = areequal(yOff,0.5*y0,1e-6,'rel');

end
