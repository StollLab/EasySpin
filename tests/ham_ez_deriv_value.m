function ok = ham_ez_deriv_value()

% Test whether magnetic moment and derivatives are correct,
% for a system with one electron

Sys.S = 1/2;
Sys.g = 2*[1 1 1];

[muMx,muMy,muMz,dmuMx,dmuMy,dmuMz] = ham_ez(Sys);


% Calculate hamiltonian from derivatives
mu = muMx + muMy + muMz;
mu_2 = Sys.g(1)*dmuMx{1}{1} + Sys.g(2)*dmuMy{1}{2} + Sys.g(3)*dmuMz{1}{3};

ok(1) = areequal(mu_2,mu,1e-10,'abs');

% test for a larger spin system and different euler angles
Sys.S = [1/2 1/2];
Sys.g = [2*[1 1.9 1.1]; 3*[1 1 2]];
Sys.gFrame = [35 26 0; 81 23 60]*pi/180;
Sys.J=1;

[muMx,muMy,muMz,dmuMx,dmuMy,dmuMz] = ham_ez(Sys);

for eSp=1:numel(Sys.S)
  ang = Sys.gFrame(eSp,:);
  if any(ang)
    R_M2g = erot(ang);  % mol frame -> g frame
  else
    R_g2M=eye(3);
  end
    R_g2M = R_M2g.';  % g frame -> mol frame
    gx{eSp} = Sys.g(eSp,1)*R_g2M(:,1)*R_g2M(:,1).';
    gy{eSp} = Sys.g(eSp,2)*R_g2M(:,2)*R_g2M(:,2).';
    gz{eSp} = Sys.g(eSp,3)*R_g2M(:,3)*R_g2M(:,3).';
end


% Calculate magnetic moment from derivatives
muMx_2 = gx{1}(1,1)*dmuMx{1}{1} + gx{1}(1,2)*dmuMx{1}{2} + gx{1}(1,3)*dmuMx{1}{3} +...
         gx{2}(1,1)*dmuMx{2}{1} + gx{2}(1,2)*dmuMx{2}{2} + gx{2}(1,3)*dmuMx{2}{3} +...
         gy{1}(1,1)*dmuMy{1}{1} + gy{1}(1,2)*dmuMy{1}{2} + gy{1}(1,3)*dmuMy{1}{3} +...
         gy{2}(1,1)*dmuMy{2}{1} + gy{2}(1,2)*dmuMy{2}{2} + gy{2}(1,3)*dmuMy{2}{3} +...
         gz{1}(1,1)*dmuMz{1}{1} + gz{1}(1,2)*dmuMz{1}{2} + gz{1}(1,3)*dmuMz{1}{3} +...
         gz{2}(1,1)*dmuMz{2}{1} + gz{2}(1,2)*dmuMz{2}{2} + gz{2}(1,3)*dmuMz{2}{3};

ok(2) = areequal(muMx_2,muMx,1e-10,'abs');

muMy_2 = gx{1}(2,1)*dmuMx{1}{1} + gx{1}(2,2)*dmuMx{1}{2} + gx{1}(2,3)*dmuMx{1}{3} +...
         gx{2}(2,1)*dmuMx{2}{1} + gx{2}(2,2)*dmuMx{2}{2} + gx{2}(2,3)*dmuMx{2}{3} +...
         gy{1}(2,1)*dmuMy{1}{1} + gy{1}(2,2)*dmuMy{1}{2} + gy{1}(2,3)*dmuMy{1}{3} +...
         gy{2}(2,1)*dmuMy{2}{1} + gy{2}(2,2)*dmuMy{2}{2} + gy{2}(2,3)*dmuMy{2}{3} +...
         gz{1}(2,1)*dmuMz{1}{1} + gz{1}(2,2)*dmuMz{1}{2} + gz{1}(2,3)*dmuMz{1}{3} +...
         gz{2}(2,1)*dmuMz{2}{1} + gz{2}(2,2)*dmuMz{2}{2} + gz{2}(2,3)*dmuMz{2}{3};

ok(3) = areequal(muMy_2,muMy,1e-10,'abs');

muMz_2 = gx{1}(3,1)*dmuMx{1}{1} + gx{1}(3,2)*dmuMx{1}{2} + gx{1}(3,3)*dmuMx{1}{3} +...
         gx{2}(3,1)*dmuMx{2}{1} + gx{2}(3,2)*dmuMx{2}{2} + gx{2}(3,3)*dmuMx{2}{3} +...
         gy{1}(3,1)*dmuMy{1}{1} + gy{1}(3,2)*dmuMy{1}{2} + gy{1}(3,3)*dmuMy{1}{3} +...
         gy{2}(3,1)*dmuMy{2}{1} + gy{2}(3,2)*dmuMy{2}{2} + gy{2}(3,3)*dmuMy{2}{3} +...
         gz{1}(3,1)*dmuMz{1}{1} + gz{1}(3,2)*dmuMz{1}{2} + gz{1}(3,3)*dmuMz{1}{3} +...
         gz{2}(3,1)*dmuMz{2}{1} + gz{2}(3,2)*dmuMz{2}{2} + gz{2}(3,3)*dmuMz{2}{3};

ok(4) = areequal(muMz_2,muMz,1e-10,'abs');


ok=all(ok);