function ok = test()

% Adding nuclei with symmetric hyperfine/quadrupole matrices [xx yy zz xy xz yz]
% Mixed forms are converted to the most general form present.

mat2sym = @(M)[M(1,1) M(2,2) M(3,3) M(1,2) M(1,3) M(2,3)];
sym2mat = @(v)[v(1) v(4) v(5); v(4) v(2) v(6); v(5) v(6) v(3)];
noFrame = @(S,F)~isfield(S,F) || isempty(S.(F));

Asym = [3 4 5 0.5 -0.2 0.1];
Apv = [1 2 6];
ang = [0.3 0.7 -0.4];
R = erot(ang).';
Apv_sym = mat2sym(R*diag(Apv)*R.');

% no prior nucleus
Sys = struct('S',1/2);
Sys = nucspinadd(Sys,'1H',Asym);
ok(1) = isequal(Sys.A,Asym);

% principal values + frame, add symmetric
Sys = struct('S',1/2,'Nucs','1H','A',Apv,'AFrame',ang);
Sys = nucspinadd(Sys,'1H',Asym);
ok(2) = areequal(Sys.A,[Apv_sym; Asym],1e-12,'abs') && noFrame(Sys,'AFrame');

% isotropic (row, two nuclei), add symmetric
Sys = struct('S',1/2,'Nucs','1H,1H','A',[1 2]);
Sys = nucspinadd(Sys,'1H',Asym);
ok(3) = areequal(Sys.A,[1 1 1 0 0 0; 2 2 2 0 0 0; Asym],1e-12,'abs');

% symmetric, add principal values + frame
Sys = struct('S',1/2,'Nucs','1H','A',Asym);
Sys = nucspinadd(Sys,'1H',Apv,ang);
ok(4) = areequal(Sys.A,[Asym; Apv_sym],1e-12,'abs') && noFrame(Sys,'AFrame');

% symmetric, add symmetric
Sys = struct('S',1/2,'Nucs','1H','A',Asym);
Sys = nucspinadd(Sys,'1H',2*Asym);
ok(5) = isequal(Sys.A,[Asym; 2*Asym]);

% symmetric, add full
Afull = magic(3);
Sys = struct('S',1/2,'Nucs','1H','A',Asym);
Sys = nucspinadd(Sys,'1H',Afull);
ok(6) = isequal(Sys.A,[sym2mat(Asym); Afull]);

% full, add symmetric
Sys = struct('S',1/2,'Nucs','1H','A',Afull);
Sys = nucspinadd(Sys,'1H',Asym);
ok(7) = isequal(Sys.A,[Afull; sym2mat(Asym)]);

% Q as [eeQq/h eta] (I=1), add symmetric Q
Qsym = [-1 -2 3 0.2 0.1 0.3];
Sys = struct('S',1/2,'Nucs','14N','A',Apv,'Q',[2 0.3]);
Sys = nucspinadd(Sys,'14N',Apv,[],Qsym);
Qpv = 2/4*[-1+0.3, -1-0.3, 2];
ok(8) = areequal(Sys.Q,[Qpv 0 0 0; Qsym],1e-12,'abs') && noFrame(Sys,'QFrame');

% symmetric matrix with nonzero frame is rejected
try
  nucspinadd(struct('S',1/2),'1H',Asym,ang);
  ok(9) = false;
catch
  ok(9) = true;
end

% 6 elements given as a 2x3 array are rejected
try
  nucspinadd(struct('S',1/2),'1H',reshape(Asym,2,3));
  ok(10) = false;
catch
  ok(10) = true;
end
