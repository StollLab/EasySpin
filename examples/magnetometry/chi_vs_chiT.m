% spin 2 with D of 1 cm-1 and axial g
Sys1 = struct('S', 2, 'g', [2.1 1.9], 'D',clight*1e-4); 
% spin 2 with D of 1 cm-1 and isotropic g
Sys2 = struct('S', 2, 'g', 1.9, 'D',clight*1e-4);       

Exp.Temperature =1:300;
Opt.Output = 'ChiTCGs ChiCGS 1overChiCGS';
[chizzT1,chizz1, OneOverChi1] = curry(Sys1,Exp,Opt);
[chizzT2,chizz2, OneOverChi2] = curry(Sys2,Exp,Opt);

xaxis = 'Temperature (K)';
leg = {'g_{\perp} = 2.1, g_{||} = 1.9', 'g_{iso} = 1.9'};
subplot(3,1,1)
plot(Exp.Temperature,[chizzT1;chizzT2])
legend(leg,'Location','best');
xlabel(xaxis); ylabel('\chi T (cm^{3} K mol^{-1})');
subplot(3,1,2)
plot(Exp.Temperature,[chizz1;chizz2])
legend(leg,'Location','best');
xlabel(xaxis); ylabel('\chi (cm^{3} mol^{-1})');
subplot(3,1,3)
plot(Exp.Temperature,[OneOverChi1;OneOverChi2])
legend(leg,'Location','best');
xlabel(xaxis); ylabel('1/\chi (cm^{-3} mol^{1})');