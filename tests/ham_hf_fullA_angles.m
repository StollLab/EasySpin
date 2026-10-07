function ok = test()

% full hyperfine matrices cannot be combined with nonzero Euler angles

Sys.S = 1/2;
Sys.Nucs = '1H,14N';
Sys.A = [rand(3); rand(3)];
Sys.AFrame = [0 0 0; rand(1,3)*pi];

try
  ham_hf(Sys);
  ok = false;
catch
  ok = true;
end
