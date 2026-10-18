function ok = test()

% Check that field modulation (Exp.ModAmp) raises an error in salt

Sys.Nucs = '1H';
Sys.A = [2 3 4];
Sys.lwEndor = 0.1;
Exp.Field = 350;
Exp.Range = [10 20];
Exp.ModAmp = 0.1;

try
  salt(Sys,Exp);
  ok = false;
catch err
  ok = contains(err.message,'salt does not support field modulation');
end

end
