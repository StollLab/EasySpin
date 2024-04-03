% mdtraj2oripot    Convert MD trajectory to orientational potential
%
%  [pot,pdf] = mdtraj2oripot(FrameTraj)
%  ___ = mdtraj2oripot(FrameTraj,nBins)
%
%  This function converts the MD frame trajectory into a orientational
%  potential energy function, defined over three Euler angles.
%
%  Input:
%     FrameTraj   MD frame trajectory
%     nBins       number-1 of bins along first and third Euler angle; second
%                 Euler angle uses half the number. For example, nBins=180
%                 gives a 181x91x181 array. Default value is 90.
%                 Alternatively, provide a vector with 3 values:
%                   [nBinsPhi, nBinsTheta, nBinsPsi]
%  Output:
%     pot         3D histogram of orientational potential energy over
%                 [phi,theta,psi]
%     pdf         3D histogram of probability distribution function over
%                 [phi,theta,psi]

function [potential,pdf] = mdtraj2oripot(FrameTraj,nBins)

if nargin<2
  nBins = 90;
end

if numel(nBins)==1
  nBins = [nBins nBins/2 nBins];
end

if numel(nBins)~=3
  error('Second input (nBins) must be a single number of an array of three numbers.');
end

% Convert frame trajectory to Euler angle trajectory
qTemp = squeeze(rotmat2quat(FrameTraj));
[phi,theta,psi] = quat2euler(qTemp,'active');
clear qTemp

% Wrap angles into [0...2*pi) interval
phi = mod(phi,2*pi);
psi = mod(psi,2*pi);

EulerTraj = [phi(:),theta(:),psi(:)];

% Build histogram over Euler angle space
phiEdges = linspace(0, 2*pi, nBins(1)+1);
thetaEdges = linspace(0, pi, nBins(2)+1);
psiEdges = linspace(0, 2*pi, nBins(3)+1);
edges = {phiEdges,thetaEdges,psiEdges};

pdf = histcountsn(EulerTraj,edges);

% Apply low-pass filter
kernel = [3 3 3];
pdf = smooth3(pdf,'gaussian',kernel);

% Put a finite floor on bin populations
% This prevents Inf and large values after taking log
thr = 1e-14;
pdf(pdf<thr) = thr;

potential = -log(pdf);

end
