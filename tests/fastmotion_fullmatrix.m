function ok = test()

% Full and symmetric A/Q matrices, and symmetric g matrices, are rejected

Sys.g = [2.0088 2.0064 2.0027];
Sys.Nucs = '14N';
Sys.A = [17 17 89];

Sys_ = Sys; Sys_.g = [Sys.g 0 0 0];
ok(1) = throwserror(Sys_);

Sys_ = Sys; Sys_.A = diag(Sys.A);
ok(2) = throwserror(Sys_);

Sys_ = Sys; Sys_.A = [Sys.A 1 0 0];
ok(3) = throwserror(Sys_);

Sys_ = Sys; Sys_.Q = diag([-1 -1 2]);
ok(4) = throwserror(Sys_);

Sys_ = Sys; Sys_.Q = [-1 -1 2 0.1 0 0];
ok(5) = throwserror(Sys_);

% full g is supported
Sys_ = Sys; Sys_.g = diag(Sys.g);
ok(6) = areequal(fastmotion(Sys_,350,1e-10),fastmotion(Sys,350,1e-10),1e-10,'rel');

end

function err = throwserror(Sys)
try
  fastmotion(Sys,350,1e-10);
  err = false;
catch
  err = true;
end
end
