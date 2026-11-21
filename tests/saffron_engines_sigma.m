function ok = test() %#ok<FNDEF>

% Chemical shielding (Sys.sigma) and gnscale in the manifold engine
% - single orientation: the manifold and full engines agree for 3p
%   ESEEM with anisotropic shielding tensors (with Euler angles) and
%   gnscale, and the shielding changes the signal
% - a full shielding tensor gives the same result as its principal values
%   and Euler angles
% - Mims ENDOR: isotropic shielding scales the nuclear Zeeman frequency, so
%   the lines are at sigma*nuI -+ A/2
% - RF pulses (full engine): with the shielding tensor diagonal in
%   the lab frame and an isotropic A, the NMR transitions are at
%   gnscale*sigma_zz*nuI -+ A/2 and the RF nutation frequency is
%   gnscale*sigma_xx*gamma*B1, so a pi pulse there gives the same ENDOR
%   effect as for a bare proton; a nominal Flip includes gnscale
% The shielding values are far larger than real ones, to make the effects
% clearly visible.

Sys.g = gfree;
Sys.Nucs = '1H,13C';
Sys.A = [0.5 1 2; -1 -0.5 1.5];  % MHz
Sys.AFrame = [0 0 0; 0.4 0.9 0.2];
Sys.sigma = [0.8 1 1.3; 1.2 0.9 1];
Sys.sigmaFrame = [0.3 0.5 0.1; 0 0.7 0.4];
Sys.gnscale = [1 1.1];

Sys0 = rmfield(Sys,{'sigma','sigmaFrame','gnscale'});

Exp.mwFreq = 9.8;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
Exp.MolFrame = [0.3 0.7 0.1];
Exp.CrystalSymmetry = 1;
Exp.Sequence = {struct('Flip',pi/2) 0.15 struct('Flip',pi/2) 0.1 struct('Flip',pi/2) 0.15};
Exp.Dim1 = {'d2' 0.05};  % µs
Exp.nPoints = 32;
Exp.DetPhase = -pi/2;

% manifold vs full engine; the full engine selects the
% stimulated echo by a phase cycle of the first and third pulses
yMan = manifold(Sys,Exp);
[ph1,ph3] = ndgrid((0:3)*pi/2,[0 pi]);
E = Exp;
E.PhaseCycle = [ph1(:) 0*ph1(:) ph3(:) -ph1(:)-ph3(:)];
[~,yFull] = saffron(Sys,E,struct('Method','full'));
ok(1) = areequal(yMan,yFull(:),1e-5,'rel');

% the shielding changes the signal much more than the tolerance
ok(2) = norm(yMan-manifold(Sys0,Exp))/norm(yMan) > 1e-3;

% full tensors in the molecular frame instead of principal values and frames
SysFull = Sys;
SysFull.sigma = zeros(6,3);
for iNuc = 1:2
  R = erot(Sys.sigmaFrame(iNuc,:));  % molecular frame -> CS frame
  SysFull.sigma(3*(iNuc-1)+(1:3),:) = R.'*diag(Sys.sigma(iNuc,:))*R;
end
SysFull = rmfield(SysFull,'sigmaFrame');
ok(3) = areequal(manifold(SysFull,Exp),yMan,1e-10,'rel');

% Mims ENDOR with isotropic shielding: lines at sigma*nuI -+ A/2
SysE.Nucs = '1H';
SysE.A = 6;  % MHz
SysE.sigma = 0.9;
SysE.lwEndor = 0.05;  % MHz
ExpE.Field = 350;  % mT
ExpE.Sequence = 'MimsENDOR';
ExpE.tau = 0.1;  % µs
ExpE.rfFreq = [6 20];  % MHz
ExpE.nPoints = 701;
[x,y] = saffron(SysE,ExpE,struct('GridSize',7));
endor = real(y(1)) - real(y);  % reduction of the echo by the RF pulse
peaks = x(islocalmax(endor) & endor>0.3*max(endor));
nuI = SysE.sigma*larmorfrq('1H',ExpE.Field);
ok(4) = numel(peaks)==2 && all(abs(peaks-(nuI+[-1 1]*SysE.A/2))<0.05);

% RF pulses with shielding and gnscale (full engine), single
% orientation with the molecular frame equal to the lab frame
SysR.Nucs = '1H';
SysR.A = 3;  % MHz
SysR.sigma = [1.3 1 0.8];  % sigma_xx, sigma_yy, sigma_zz
SysR.gnscale = 1.1;
SysRef.Nucs = '1H';
SysRef.A = SysR.A;
gamma = nucgval('1H')*nmagn/planck/1e9;  % MHz/mT
p90.Flip = pi/2;
rf.Channel = 'rf';
rf.tp = 2;  % µs
ExpR.Field = 350;  % mT
ExpR.mwFreq = 9.8;  % GHz
ExpR.MolFrame = [0 0 0];
ExpR.CrystalSymmetry = 1;
ExpR.Pathways = '+00-';
ExpR.DetPhase = -pi/2;
nuLine = @(s) s*gamma*ExpR.Field + SysR.A/2;  % upper NMR transition, MHz
B1pi = @(s) 1/(2*rf.tp*s*gamma);  % RF field (mT) for a pi pulse
eRef = rfEffect(SysRef,ExpR,rf,p90,nuLine(1),B1pi(1));
eR = rfEffect(SysR,ExpR,rf,p90,nuLine(1.1*0.8),B1pi(1.1*1.3));
ok(5) = abs(eRef)>0.1 && areequal(eR,eRef,1e-3,'rel');

% nominal Flip = pi refers to the first nucleus, with gnscale
SysG = SysRef;
SysG.gnscale = 1.1;
eG = rfEffect(SysG,ExpR,rf,p90,nuLine(1.1),[]);
ok(6) = areequal(eG,eRef,1e-3,'rel');

end

% ENDOR effect of an RF pulse at frequency nuRF (MHz), with RF field B1
% (mT) or, for empty B1, Flip = pi: echo without minus echo with RF field,
% in Mims ENDOR (full engine)
function e = rfEffect(Sys,Exp,rf,p90,nuRF,B1)
if isempty(B1)
  rf.Flip = pi;
  rfOff = rf; rfOff.Flip = 0;
else
  rf.B1 = B1;
  rfOff = rf; rfOff.B1 = 0;
end
Exp.rfFreq = nuRF;
Opt.Method = 'full';
Exp.Sequence = {p90 0.2 p90 0 rf 0.1 p90 0.2};
[~,yOn] = saffron(Sys,Exp,Opt);
Exp.Sequence{5} = rfOff;
[~,yOff] = saffron(Sys,Exp,Opt);
e = real(yOff-yOn);
end

function y = manifold(Sys,Exp)
Exp.Pathways = 'echo';
[~,y] = saffron(Sys,Exp,struct('Method','manifold'));
y = y(:);
end
