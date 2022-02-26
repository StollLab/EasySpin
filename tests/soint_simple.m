function ok = test()

% Test against explicit values for simple case

Sys.S = 1/2;
Sys.L = 1;
Sys.soc = 1;

H = soint(Sys);

a = 0.5;
b = 1/sqrt(2);
Href = ...
  [a 0 0 0 0 0;
    0 0 0 b 0 0;
    0 0 -a 0 b 0;
    0 b 0 -a 0 0;
    0 0 b 0 0 0;
    0 0 0 0 0 a];

ok = areequal(H,Href,1e-10,'abs');
