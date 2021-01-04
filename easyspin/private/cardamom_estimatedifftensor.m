%  cardamom_estimatedifftensor Estimate the rotational diffusion tensor
%                              from an MD simulation trajectory.
%
%  [Diff,msadp,tLag] = cardamom_estimatedifftensor(RTraj_L, t);
%
%  Input:
%      RTraj          numeric, size = (3,3,nSteps)
%                     rotation matrices in the lab frame
%      dt             double
%                     time step
%      stopFitT       double
%                     time value at which to stop the fit, depending on the
%                     timescales of the MD and EPR simulations
%
%  Output:
%      Diff           numeric, size = (3,3)
%                     rotational diffusion tensor
%      msadp          numeric, size = (3,round(end/2))
%                     mean square angular displacement in the body-fixed
%                     frame
%      tLag           double
%                     time lag

% Implemented from
%    [1] G. Chevrot, et al., J. Chem. Phys. 139, 154110 (2013)
%        http://dx.doi.org/10.1063/1.4823996
%    [2] V. Calandrini, et al., Collection SFN 12, 201 (2011)
%        	https://doi.org/10.1051/sfn/201112010


function [Diff, msadp, tLag] = cardamom_estimatedifftensor(RTraj, dt, stopFitT)

nSteps = length(RTraj);

% orient the trajectory such that the tensor of inertia is diagonal at the
% first time point (part of Step 1 in Sec. IIIB of Ref. 1)
RAlign = RTraj(:,:,1);
RTrajp = zeros(3,3,nSteps);
for iStep = 1:nSteps
  RTrajp(:,:,iStep) = RTraj(:,:,iStep)*RAlign.';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate rotation matrices that transform (propagate) the molecular
% coordinates between time steps
RRot = zeros(3,3,nSteps-1);
for iStep = 1:nSteps-1
  RRot(:,:,iStep) = RTrajp(:,:,iStep+1)*RTrajp(:,:,iStep).';
end

% convert to quaternions
qRot = rotmat2quat(RRot);
qRot = cat(2, [1;0;0;0], qRot);  % first time point is the identity

% ensure that the quaternion first component, which is associated with the
% angular part of the axis-angle formulation, is always positive, thereby
% restricting to one half of the hypersphere
idx = qRot(1,:) < 0;
qRot(:,idx) = -qRot(:,idx);

% calculate Cartesian angular velocity components in the body-fixed frame
wp = q2wp(qRot, dt);

% calculate angular displacement in the body-fixed frame
t = dt*(0:length(wp)-1);
Deltawp = cumtrapz(t, wp, 2);  % Eq. 32 of Ref. 1

% calculate mean square angular displacement
msadp = msd_fft(Deltawp);
msadp = msadp(:, 1:round(end/2));  % statistics for the latter half of the
                                   % MSAD are poor, so only keep the first
                                   % half


% estimate the eigenvalues of the rotational diffusion tensor using
% least-squares fitting
tLag = t(1:length(msadp))/1e-9;
stopFitN = floor(stopFitT/dt);
pxp = polyfit(tLag(1:stopFitN), msadp(1,1:stopFitN), 1);
pyp = polyfit(tLag(1:stopFitN), msadp(2,1:stopFitN), 1);
pzp = polyfit(tLag(1:stopFitN), msadp(3,1:stopFitN), 1);

Diff = [pxp(1), pyp(1), pzp(1)]/2*1e9;

end

% Helper functions
% -------------------------------------------------------------------------

function w = q2w(qTraj, dt)
% Convert quaternion trajectory to angular velocity in the lab frame
% See Eq. 29 of Ref. 1

dq = diff(qTraj, 1, 2)/dt;

q0 = qTraj(1, 1:end-1);
q1 = qTraj(2, 1:end-1);
q2 = qTraj(3, 1:end-1);
q3 = qTraj(4, 1:end-1);

dq0 = dq(1, :);
dq1 = dq(2, :);
dq2 = dq(3, :);
dq3 = dq(4, :);

wx = 2*(-q1.*dq0 + q0.*dq1 - q3.*dq2 + q2.*dq3);
wy = 2*(-q2.*dq0 + q3.*dq1 + q0.*dq2 - q1.*dq3);
wz = 2*(-q3.*dq0 - q2.*dq1 + q1.*dq2 + q0.*dq3);

w = [wx; wy; wz];

end

function wp = q2wp(qTraj, dt)
% Convert quaternion trajectory to angular velocity in the molecule-fixed
% frame
% See Eq. 30 of Ref. 1

dq = diff(qTraj, 1, 2)/dt;

q0 = qTraj(1, 1:end-1);
q1 = qTraj(2, 1:end-1);
q2 = qTraj(3, 1:end-1);
q3 = qTraj(4, 1:end-1);

dq0 = dq(1, :);
dq1 = dq(2, :);
dq2 = dq(3, :);
dq3 = dq(4, :);

wxp = 2*(-q1.*dq0 + q0.*dq1 + q3.*dq2 - q2.*dq3);
wyp = 2*(-q2.*dq0 - q3.*dq1 + q0.*dq2 + q1.*dq3);
wzp = 2*(-q3.*dq0 + q2.*dq1 - q1.*dq2 + q0.*dq3);

wp = [wxp; wyp; wzp];

end

function msd = msd_fft(x)
% calculate the mean square displacement using the FFT
% see Sec. 4.2 in Ref. 2

if iscolumn(x)
  x = x.';
end

nComps = size(x, 1);
N = length(x);

D = zeros(nComps, N+1);
D(:,2:end) = x.^2;


S2 = runprivate('autocorrfft', x, 2, 0, 0, 0);

Q = 2*sum(D, 2);
S1 = zeros(nComps, N);

for m = 1:N
    Q = Q - D(:, m) - D(:, end-m);
    S1(:, m) = Q/((N+1)-m);
end

msd = S1 - 2*S2;

end
