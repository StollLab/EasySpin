% orifun_M2L_chiint  Evaluates ordering function
%
%  w = orifun_M2L_chiint(phi,theta,R_L2S,orifun_S2M)
%
% Inputs:
%   phi    phi angle for the transformation from mol to lab frame
%   theta  theta angle for the transformation from mol to lab frame
%   R_L2S  transformation matrix from lab to sample frame
%   f      orientational distribution function, taking as input
%          Euler angles from sample frame to molecular frame
%
% Output:
%   w      weights, one for each element in phi/theta,
%          integrated over chi

function w = orifun_M2L_chiint(phi,theta,R_L2S,f)

% Set up grid over third angle, chi
nChi = 12;
chi = linspace(0,2*pi,nChi);
dchi = chi(2)-chi(1);
c = cos(chi);
s = sin(chi);

R_S2L = R_L2S.';
for iOri = numel(phi):-1:1
  R_M2L = erot(phi(iOri),theta(iOri),0);
  R_L2M = R_M2L.';
  xL_M = R_L2M(:,1)*c + R_L2M(:,2)*s;
  yL_M = -R_L2M(:,1)*s + R_L2M(:,2)*c;
  for iChi = numel(chi):-1:1
    R_L2M(:,1) = xL_M(:,iChi);
    R_L2M(:,2) = yL_M(:,iChi);
    R_S2M = R_L2M*R_S2L;
    [alpha(iChi),beta(iChi),gamma(iChi)] = eulang(R_S2M,true);
  end
  ow = f(alpha,beta,gamma);
  w(iOri) = trapz(ow)*dchi;
end

end
