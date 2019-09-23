function [err,data] = test(opt,olddata)

% Test whether salt() can handle systems
% with a mixture of magnetic and non-magnetic isotopes

Sys.Nucs = 'Ti';
Sys.A = 10;
Sys.lwEndor = 1;

Exp.Field = 350;
Exp.Range = [0 50];

[rf,spec] = salt(Sys,Exp);

if opt.Display
  plot(rf,spec,rf,olddata.spec);
  legend('new','old')
end

if isempty(olddata)
  err = [];
else
  err = ~areequal(spec,olddata.spec,1e-10,'abs');
end

data.spec = spec;
