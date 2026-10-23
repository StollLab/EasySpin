function ok = test()

% Hyperfine enhancement is rejected with perturbation theory

Sys = struct('S',1/2,'g',[2 2.1 2.2],'Nucs','1H','A',[3 7 12]);
Exp = struct('Field',350,'SampleFrame',[pi/5 pi/7 pi/6]);

Opt.Enhancement = 'off';
endorfrq_perturb(Sys,Exp,Opt);

Opt.Enhancement = 'on';
try
  endorfrq_perturb(Sys,Exp,Opt);
  ok = false;
catch err
  ok = contains(err.message,'Hyperfine enhancement');
end
