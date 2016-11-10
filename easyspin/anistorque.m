%  anistorque  Calculate the torque due to an anisotripic orienting potential.
%
%   iJdotu = anistorque(Q, coef);
%
%   Input:
%      Q              2x2xnxm array representation of a quaternion
%      coefs          a struct containing the orienting potential
%                     coefficients, field values should be specified in the
%                     format 'cjpm', 'cjmm', or 'cj0'
%                     
%
%   Output:
%      anistorque     3x1xn vector

% Implementation based on 
% - Sezer, et al., J.Chem.Phys. 128, 165106 (2008), doi: 10.1063/1.2908075

function torque = anistorque(Q, coefs)

sizeQ = size(Q);
shapeQ = num2cell(sizeQ);
Index = cell(1, ndims(Q));
Index(:) = {':'};

if ~isnumeric(Q) || sizeQ(1)~=2 || sizeQ(2)~=2
  error('The first two dimensions of the array must be of shape 2x2.')
end

if ~isstruct(coefs)
  error('Orienting potential coefficients need to be stored in a structured array.')
end

fields = fieldnames(coefs);
% Need to check formatting of fields!

% jVals = cellfun(@(x) str2double(x(2)), fields);
% mVals = cellfun(@(x) str2double(x(3)), fields);

Nfields = numel(fields);
m1 = 0;
if ndims(Q)>2
  iJvecu = zeros(3,shapeQ{3:end},Nfields);
else
  iJvecu = zeros(3,Nfields);
end

% Not sure if there is a good way to vectorize this
for ifield=1:Nfields
  fieldstr = fields{ifield};
  j = str2double(fieldstr(2));
  if fieldstr(3)=='p'
    m2 = str2double(fieldstr(4));
  elseif fieldstr(3)=='m'
    m2 = -str2double(fieldstr(4));
  elseif fieldstr(3)=='0'
    m2=0;
  else
    error('Potential coefficient formatting is not correct.')
  end
  coef = coefs.(fieldstr);
  iJvecu(:,Index{3:end},ifield) = iJu(j,m1,m2,coef,Q);
end

if Nfields > 1
  torque = -sum(iJvecu,ndims(iJvecu));
else
  torque = -iJvecu;
end

% torque = -sum(iJvecu,2);

% Q11 = Q(1,1,:);
% Q12 = Q(1,2,:);
% Q21 = Q(2,1,:);
% Q22 = Q(2,2,:);

% % Eq. C7 and C8 in reference
% p_ = Q11.*Q22 + Q12.*Q21;
% D2_0p1 = sqrt(6)*Q11.*Q21.*p_;
% D2_0m1 = sqrt(6)*Q22.*Q12.*p_;
% 
% iJxu = 1i*sqrt(3/2)*c20*(D2_0p1 + D2_0m1);
% iJyu =    sqrt(3/2)*c20*(D2_0p1 - D2_0m1);
% iJzu = zeros(size(iJyu));
% 
% torque = -[iJxu; iJyu; iJzu];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  out = iJu(j, m1, m2, cjm, Q)
% Calculate the result of i\vo{J}u for a given potential coefficient
% cjm corresponds to c_{m}^{j} in Eq. 58 in reference

if m2==j
  % m2 cannot be greater than j, so D_{0,m2+1}^j = 0
  iJxu = 1i/2*cjm*Cpm(j,m2,'-')*wignerdquat(j,m1,m2-1,Q);
  iJyu = -1/2*cjm*Cpm(j,m2,'-')*wignerdquat(j,m1,m2-1,Q);
elseif m2==-j
  % m2 cannot be less than -j, so D_{0,m2-1}^j = 0
  iJxu = 1i/2*cjm*Cpm(j,m2,'+')*wignerdquat(j,m1,m2+1,Q);
  iJyu = 1/2*cjm*Cpm(j,m2,'+')*wignerdquat(j,m1,m2+1,Q);
else
  iJxu = 1i/2*cjm*(Cpm(j,m2,'+')*wignerdquat(j,m1,m2+1,Q) ...
                  +Cpm(j,m2,'-')*wignerdquat(j,m1,m2-1,Q));
  iJyu = 1/2*cjm*(Cpm(j,m2,'+')*wignerdquat(j,m1,m2+1,Q) ...
                 -Cpm(j,m2,'-')*wignerdquat(j,m1,m2-1,Q));
end

iJzu = 1i*cjm*m2*wignerdquat(j,m1,m2,Q);    

out = [iJxu; iJyu; iJzu];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = Cpm(j, m, pm)
% Eq. 60 in reference

if pm=='+', C = sqrt(j*(j+1)-m*(m+1));
elseif pm=='-', C = sqrt(j*(j+1)-m*(m-1)); end

end

end
