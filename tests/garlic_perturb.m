function ok = test(opt)

% Test various computation methods (perturbation, fixed-point)

Sys = struct('g',2,'Nucs','1H','A',600,'n',3,'lw',[0 0.1]);
Exp = struct('mwFreq',9.7,'nPoints',10000,'Range',[310 380]);
Sim.Method = 'perturb1';
[x,y(1,:)] = garlic(Sys,Exp,Sim);
Sim.Method = 'perturb2';
[x,y(2,:)] = garlic(Sys,Exp,Sim);
Sim.Method = 'perturb3';
[x,y(3,:)] = garlic(Sys,Exp,Sim);
Sim.Method = 'perturb4';
[x,y(4,:)] = garlic(Sys,Exp,Sim);
Sim.Method = 'perturb5';
[x,y(5,:)] = garlic(Sys,Exp,Sim);
Sim.Method = 'exact';
[x,y(6,:)] = garlic(Sys,Exp,Sim);

if opt.Display
  plot(x,y); xlabel('magnetic field [mT]');
  axis tight, set(gca,'YTick',[]);
  legend('1st order','2nd order','3rd order','4th order','5th order','Breit-Rabi');
  pause;
end

ok = true;
