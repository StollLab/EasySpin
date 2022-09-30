% mdtraj2oripot    Convert MD trajectory to orientational potential
%
%  [pot,pdf] = mdtraj2oripot(FrameTraj)
%  ___ = mdtraj2oripot(FrameTraj,nBins)
%
%  This function converts the MD frame trajectory into a orientational
%  potential energy function, defined over three Euler angles
%
%  Input:
%     FrameTraj   MD frame trajectory
%     nBins       number-1 of bins along first and third Euler angle; second
%                 Euler angle uses half the number. For example, nBins=180
%                 gives a 181x91x181 array. Default value is 90.
%  Output:
%     pot         3D histogram of orientational potential energy
%     pdf         3D histogram of probability distribution function

function [potential,pdf] = mdtraj2oripot(FrameTraj,nBins)

if nargin<2, nBins = 90; end

% Convert frame trajectory to Euler angle trajectory
qTemp = squeeze(rotmat2quat(FrameTraj));
[phi,theta,psi] = quat2euler(qTemp,'active');

% Replace negative Euler angles with equivalent positive ones
phi = phi + 2*pi*(phi<0);
psi = psi + 2*pi*(psi<0);

EulerTraj = [phi(:),theta(:),psi(:)];
clear qTemp

% Build histogram over Euler angle space
phiEdges = linspace(0, 2*pi, nBins+1);
thetaEdges = linspace(0, pi, nBins/2+1);
psiEdges = linspace(0, 2*pi, nBins+1);
edges = {phiEdges,thetaEdges,psiEdges};
pdf = histcountsn(EulerTraj,edges);

% Apply low-pass filter
pdf = smooth3(pdf,'gaussian');

% Put a finite floor on histogram (prevents inf and large values after log)
thr = 1e-14;
pdf(pdf<thr) = thr;

potential = -log(pdf);

end
