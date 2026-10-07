function ok = test()

% full g matrices cannot be combined with nonzero Euler angles

B = rand(1,3)*340;

Sys.S = 3/2;
Sys.Nucs = '1H';
Sys.A = [30 40 50];
Sys.g = rand(3);
Sys.gFrame = rand(1,3)*2*pi;

try
  ham_ez(Sys,B);
  ok = false;
catch
  ok = true;
end
