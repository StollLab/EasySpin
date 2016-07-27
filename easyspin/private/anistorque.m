%  anistorque  Calculate the torque due to an anisotripic orienting potential.
%
%  iJdotu = anistorque(Q, coef);
%
%  Input:
%     Q              2x2xnxm array representation of a quaternion
%     c20            float, orienting potential coeficient
%                    Only c20 is currently implemented!
%
%  Output:
%     anistorque     3x1xn vector

% Implementation based on 
% - Sezer, et al., J.Chem.Phys. 128, 165106 (2008), doi: 10.1063/1.2908075

function torque = anistorque(Q, c20)

if size(Q,1)~=2 || size(Q,2)~=2
  error('The first two dimensions of the array must be of shape 2x2.')
end

Q11 = Q(1,1,:);
Q12 = Q(1,2,:);
Q21 = Q(2,1,:);
Q22 = Q(2,2,:);

% Eq. C7 and C8 in reference
p_ = Q11.*Q22 + Q12.*Q21;
D2_0p1 = sqrt(6)*Q11.*Q21.*p_;
D2_0m1 = sqrt(6)*Q22.*Q12.*p_;

iJxu = 1i*sqrt(3/2)*c20*(D2_0p1 + D2_0m1);
iJyu =    sqrt(3/2)*c20*(D2_0p1 - D2_0m1);
iJzu = zeros(size(iJyu));

torque = -[iJxu; iJyu; iJzu];

end
