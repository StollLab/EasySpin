function ok = test()

% full D matrices cannot be combined with nonzero Euler angles

Sys.S = 3/2;
Sys.D = rand(3);
Sys.DFrame = rand(1,3)*2*pi;

try
  ham_zf(Sys);
  ok = false;
catch
  ok = true;
end
