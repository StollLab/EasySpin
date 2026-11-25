% Mims ENDOR of a strongly coupled 1H: hyperfine enhancement (saffron)
%==========================================================================
% The RF field acts on a nucleus through its nuclear Zeeman interaction,
% but it also tilts the electron spin slightly. Through the hyperfine
% coupling, this changes the field the nucleus feels (hyperfine
% enhancement). The RF pulse is then stronger for the nuclear transitions
% in one electron spin manifold and weaker in the other, which gives the
% two ENDOR lines of a nucleus different intensities. The effect grows with
% the ratio of the hyperfine coupling to the nuclear Zeeman frequency.
%
% By default (Opt.rfEnhancement = true), saffron includes it. Here, the
% coupling (10-16 MHz) is close to the nuclear Zeeman frequency (14.9 MHz
% at 350 mT). The RF pulse is weak (nominal flip angle pi/4 for a bare
% proton), so the transfer grows with the square of the RF field at the
% nucleus, and the enhancement raises one line and lowers the other.
%
% Reference: A. Schweiger, G. Jeschke, Principles of Pulse Electron
% Paramagnetic Resonance, Oxford University Press, 2001, Section 6.3.

clear, clc

% Spin system: S=1/2 with one strongly coupled 1H
Sys.Nucs = '1H';
Sys.A = [10 10 16];  % hyperfine coupling, MHz
Sys.lwEndor = 0.1;   % ENDOR line width, MHz

% Experiment: Mims ENDOR with a weak RF pulse
p90.Flip = pi/2;
RF.Channel = 'rf';
RF.tp = 10;          % µs
RF.Flip = pi/4;
tau = 0.12;          % µs
Exp.Field = 350;     % mT
Exp.Sequence = {p90 tau p90 0.5 RF 0.5 p90 tau};
Exp.Pathways = '+00-';  % stimulated echo
Exp.DetPhase = -pi/2;   % stimulated echo along the real axis
Exp.rfFreq = 4;      % MHz
Exp.nPoints = 441;
Exp.Dim1 = {'rfFreq',0.05};

Opt.GridSize = 91;

% Simulate with and without hyperfine enhancement
Opt.rfEnhancement = true;
[rf,yEnh] = saffron(Sys,Exp,Opt);
Opt.rfEnhancement = false;
[~,yBare] = saffron(Sys,Exp,Opt);

% ENDOR effect: change of the echo by the RF pulse
endorEnh = real(yEnh(1)-yEnh);
endorBare = real(yBare(1)-yBare);

plot(rf,endorEnh,rf,endorBare);
axis tight
grid on
xlabel('RF frequency (MHz)');
ylabel('ENDOR effect');
legend('with hyperfine enhancement','nuclear Zeeman only','Location','north');
title('Mims ENDOR, 1H with A close to \nu_I');
