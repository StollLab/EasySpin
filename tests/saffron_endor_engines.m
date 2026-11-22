function ok = test() %#ok<FNDEF>

% ENDOR with one rectangular RF pulse: the manifold engine (two-level
% population transfer per NMR transition) and the full engine
% (propagation of the RF pulse) agree, with and without the hyperfine
% enhancement (Opt.rfEnhancement), for RF pulses given by B2 or by Flip.
% Tolerances, relative to the largest ENDOR effect:
% - field along the hyperfine axis, 1H: 1e-3, including the side lobes of
%   the RF excitation profile
% - tilted hyperfine tensor, crystal and small powder: 1% and 3%. The
%   manifold engine neglects the second-order hyperfine shifts of the NMR
%   lines, and the change of nuclear coherences (from the microwave
%   pulses) by the RF pulse
% - Davies ENDOR with selective microwave pulses, tilted hyperfine tensor:
%   4%, for the same reasons, with the larger hyperfine coupling (about
%   3% here; 4e-4 with the hyperfine axis along the field)
% - 14N with A close to 2*nuI (strong enhancement), field along the tensor
%   axes: 5%, from the Bloch-Siegert shift of the lines at low frequencies
%   (rotating-wave approximation in the manifold engine)

Sys.Nucs = '1H';
Sys.A = [1 1 3];  % MHz
Sys.g = 2.0023;

p90.Flip = pi/2;
rf.Channel = 'rf';
rf.tp = 2;  % µs
rf.B2 = 0.5/(2*0.04258);  % mT, pi pulse for a bare proton

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0 0];
Exp.CrystalSymmetry = 1;
Exp.rfFreq = 13;  % MHz
Exp.Sequence = {p90 0.2 p90 0 rf 0.1 p90 0.2};
Exp.nPoints = 21;
Exp.Dim1 = {'rfFreq',0.2};
Exp.Pathways = '+00-';
Exp.DetPhase = -pi/2;

[ok(1),x,y1,y1off,y2off,info1,info2] = compareEngines(Sys,Exp,struct(),1e-3);
ok(2) = strcmp(info1.Method,'manifold') && strcmp(info2.Method,'full');
ok(3) = areequal(x,13+0.2*(0:20));
ok(4) = areequal(y1off,y2off,1e-4,'rel') && max(abs(y1off-y1))>0.1*abs(y1off);

% the sign of the RF field does not change the population transfer
Eneg = Exp;
Eneg.Sequence{5}.B2 = -rf.B2;
[~,y1neg] = saffron(Sys,Eneg,struct('Method','manifold'));
ok(5) = areequal(y1neg,y1,1e-10,'rel');

% Tilted hyperfine tensor, RF pulse with B2 and with Flip, crystal and
% small powder
SysT = Sys;
SysT.A = [-1 -1 5];  % MHz
SysT.AFrame = [0 0.4 0];
ExpT = Exp;
ExpT.MolFrame = [0 0.7 0];
ExpT.rfFreq = 11;  % MHz
ExpT.nPoints = 41;
ok(6) = compareEngines(SysT,ExpT,struct(),0.01);
ExpT.Sequence{5} = struct('Channel','rf','tp',2,'Flip',pi);
ok(7) = compareEngines(SysT,ExpT,struct(),0.01);
ExpP = rmfield(ExpT,{'MolFrame','CrystalSymmetry'});
ExpP.nPoints = 21;
ExpP.Dim1 = {'rfFreq',0.4};
ok(8) = compareEngines(SysT,ExpP,struct('GridSize',7),0.03);

% Davies ENDOR, selective microwave pulses
pS90 = struct('Flip',pi/2,'tp',0.05);
pS180 = struct('Flip',pi,'tp',0.1);
SysD = SysT;
SysD.A = [-2 -2 10];  % MHz
ExpD = ExpT;
ExpD.Sequence = {pS180 0 ExpT.Sequence{5} 0.5 pS90 0.3 pS180 0.3};
ExpD.Pathways = '00+-';
ExpD.mwFreq = 9.8109;  % GHz, 1.5 MHz above the center of the EPR line
ok(9) = compareEngines(SysD,ExpD,struct(),0.04);

% 14N, A close to 2*nuI, field along the tensor axes, lines above 3 MHz
SysN.Nucs = '14N';
SysN.A = [2 2 3];  % MHz
SysN.Q = [-1 -1 2]*0.8;  % MHz
ExpN = ExpT;
ExpN.MolFrame = [0 0 0];
ExpN.Sequence{5}.tp = 6;  % µs
ExpN.rfFreq = 3;  % MHz
ExpN.Dim1 = {'rfFreq',0.05};
ExpN.nPoints = 51;
ok(10) = compareEngines(SysN,ExpN,struct(),0.05);

end

% ENDOR effects (reduction of the echo by the RF pulse) of both engines,
% with and without hyperfine enhancement, agree within tol relative to the
% largest ENDOR effect. The echo without RF field (B2 or Flip zero) is the
% same at all RF frequencies, so a single point suffices.
function [ok,x,y1,y1off,y2off,info1,info2] = compareEngines(Sys,Exp,Opt,tol)
iRF = find(cellfun(@(p) isstruct(p) && isfield(p,'Channel'),Exp.Sequence));
Eoff = Exp;
if isfield(Eoff.Sequence{iRF},'B2')
  Eoff.Sequence{iRF}.B2 = 0;
else
  Eoff.Sequence{iRF}.Flip = 0;
end
Eoff.nPoints = 1;
ok = true;
for enh = [false true]
  Opt.rfEnhancement = enh;
  Opt.Method = 'manifold';
  [x,y1,info1] = saffron(Sys,Exp,Opt);
  [~,y1off] = saffron(Sys,Eoff,Opt);
  Opt.Method = 'full';
  [~,y2,info2] = saffron(Sys,Exp,Opt);
  [~,y2off] = saffron(Sys,Eoff,Opt);
  endor1 = real(y1off-y1);
  endor2 = real(y2off-y2);
  ok = ok && areequal(endor1,endor2,tol*max(abs(endor2)),'abs');
end
end
