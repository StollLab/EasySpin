% Speed-up of ENDOR simulations using orientation pre-selection
%==========================================================================
% Opt.OriPreSelect provides a speed-up for simulations of
% strongly orientation-selective ENDOR spectra. If the orientation
% selection is not strong, the speed-up is only minor.

clear, clf

% Spin system
Sys.g = [2,2,2.2];
Sys = nucspinadd(Sys,'63Cu',[40 40 400]);
Sys = nucspinadd(Sys,'1H',[2,3,10]);
Sys.HStrain = [1 1 1]*80;     % MHz
Sys.lwEndor = 0.1;    % MHz

% ENDOR simulations
B0 = 315; % mT
EndorExp.Field = B0;
EndorExp.Range = larmorfrq('1H',B0) + [-1 1]*10; % MHz
EndorExp.mwFreq = 9.5;
EndorExp.ExciteWidth = 50;

Opt.nKnots = 41;
Opt.Nuclei = 2;   % simulate only ENDOR of second nucleus, i.e., 1H

% (1) no orientation pre-selection: slow
Opt.OriPreSelect = 0;
tic
[x,y1] = salt(Sys,EndorExp,Opt);
toc

% (2) with orientation pre-selection: faster
Opt.OriPreSelect = 1;
tic
[x,y2] = salt(Sys,EndorExp,Opt);
toc

plot(x,y1,'g',x,y2,'r')

xlabel('frequency [MHz]');
title('1H ENDOR spectrum');
legend('no pre-selection','with pre-selection');
