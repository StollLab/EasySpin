% orifun_M2L  Evaluates ordering function
%
%  w = orifun_M2L(f,R_L2S,phi,theta)
%  w = orifun_M2L(f,R_L2S,phi,theta,chi)
%
% Inputs:
%   phi    phi angle for the transformation from mol to lab frame
%   theta  theta angle for the transformation from mol to lab frame
%   chi    chi angle for the transformation from mol to lab frame
%          if omitted, f is integrated over chi
%   R_L2S  transformation matrix from lab to sample frame
%   f      orientational distribution function, taking as input
%          Euler angles from sample frame to molecular frame
%
% Output:
%   w      weights, one for each element in phi/theta,
%          integrated over chi if chi is not given

function w = orifun_M2L(f,R_L2S,phi,theta,chi)

R_S2L = R_L2S.';

integrateChi = nargin==4;

if integrateChi

  bruteForce = false;
  if bruteForce

    for iOri = numel(phi):-1:1
      Fchi = @(chi) oriweight(f,R_S2L,phi(iOri),theta(iOri),chi);
      w(iOri) = integral(Fchi,0,2*pi);
    end

  else

    % Set up simple coarse integration grid over chi
    nChi = 12;
    chi = linspace(0,2*pi,nChi);
    dchi = chi(2)-chi(1);
    c = cos(chi);
    s = sin(chi);

    % Run over all orientations
    for iOri = numel(phi):-1:1
      R_M2L = erot(phi(iOri),theta(iOri),0);
      R_L2M = R_M2L.';
      xL_M = R_L2M(:,1)*c + R_L2M(:,2)*s;
      yL_M = -R_L2M(:,1)*s + R_L2M(:,2)*c;
      for iChi = numel(chi):-1:1
        R_L2M(:,1) = xL_M(:,iChi);
        R_L2M(:,2) = yL_M(:,iChi);
        R_S2M = R_L2M*R_S2L;
        % Calculate Euler angles from S to M frame
        [alpha_S2M(iChi),beta_S2M(iChi),gamma_S2M(iChi)] = eulang(R_S2M,true);
      end
      ow = f(alpha_S2M,beta_S2M,gamma_S2M);
      w(iOri) = trapz(ow)*dchi;
    end

  end

else

  for iOri = numel(phi):-1:1
    R_M2L = erot(phi(iOri),theta(iOri),chi(iOri));
    R_L2M = R_M2L.';
    R_S2M = R_L2M*R_S2L;
    [alpha_S2M,beta_S2M,gamma_S2M] = eulang(R_S2M,true);
    w(iOri) = f(alpha_S2M,beta_S2M,gamma_S2M);
  end

end

end

function ow = oriweight(f,R_S2L,ph,th,chi)
for iChi = numel(chi):-1:1
  R_M2L = erot(ph,th,chi(iChi));
  R_L2M = R_M2L.';
  R_S2M = R_L2M*R_S2L;
  [alpha_S2M,beta_S2M,gamma_S2M] = eulang(R_S2M,true);
  ow(iChi) = f(alpha_S2M,beta_S2M,gamma_S2M);
end
end
