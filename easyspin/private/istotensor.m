% istotensor takes two vector spin operators as input and gives as output
% the 1 isotropic and 5 rank-2 Cartesian irreducible spherical tensor
% components. The isotropic component is stored as a row of 1 matrix in a
% cell array. The anisotropic components are arranged as a row of 5
% matrices in a cell array in decreasing order of component M from 2 to -2.

function [T0,T1,T2] = istotensornew(A,B)


Ax = A{1}; Ay = A{2}; Az = A{3};
Bx = B{1}; By = B{2}; Bz = B{3};

CartesianTensors = false;

if CartesianTensors
  
  % L = 0
  %-----------------------------------------
  T0 = -(1/sqrt(3)) * (Ax*Bx+Ay*By+Az*Bz); % M = 0
  
  % L = 1
  %-----------------------------------------
  T1 = cell(1,3);
  T1{1} = -0.5*((Ax*Bz - Az*Bx) + 1i*(Ay*Bz - Az*By)); % M = +1
  T1{3} = -0.5*((Ax*Bz - Az*Bx) - 1i*(Ay*Bz - Az*By)); % M = -1
  T1{2} = -(1/sqrt(2))*(Ay*Bx - Ax*By);               % M =  0
  
  % L = 2
  %-----------------------------------------
  T2 = cell(1,5);
  T2{1} =  0.5*((Ax*Bx - Ay*By) + 1i*(Ax*By + Ay*Bx)); % M = +2
  T2{5} =  0.5*((Ax*Bx - Ay*By) - 1i*(Ax*By + Ay*Bx)); % M = -2
  T2{2} = -0.5*((Ax*Bz + Az*Bx) + 1i*(Ay*Bz + Az*By)); % M = +1
  T2{4} =  0.5*((Ax*Bz + Az*Bx) - 1i*(Ay*Bz + Az*By)); % M = -1
  T2{3} =  sqrt(2/3) *  (Az*Bz - 0.5*(Ax*Bx + Ay*By)); % M =  0
  
else  % polar representation
  
  Ar   = Ax + 1i*Ay;  % raising operators  (up)
  Br   = Bx + 1i*By;
  
  Al = Ax - 1i*Ay;    % lowering operators (down)
  Bl = Bx - 1i*By;
  
  % L = 0
  %-----------------------------------------
  T0 = -(1/sqrt(3)) * (Az*Bz + 0.5*(Ar*Bl + Al*Br)); % M = 0
  
    % L = 1
  %-----------------------------------------
  T1 = cell(1,3);
  T1{1} = -0.5 * (Ar*Bz - Az*Br);            % M = +1
  T1{3} = -0.5 * (Al*Bz - Az*Bl);            % M = -1
  T1{2} = -(0.5i/sqrt(2)) * (Ar*Bl - Al*Br); % M =  0      
  
  % L = 2
  %-----------------------------------------
  T2 = cell(1,5);
  T2{1} =  0.5 * (Ar*Br);                              % M = +2
  T2{5} =  0.5 * (Al*Bl);                              % M = -2
  T2{2} = -0.5 * (Ar*Bz + Az*Br);                      % M = +1
  T2{4} =  0.5 * (Al*Bz + Az*Bl);                      % M = -1
  T2{3} =  sqrt(2/3) * (Az*Bz - 0.25*(Ar*Bl + Al*Br)); % M =  0
  
end

return