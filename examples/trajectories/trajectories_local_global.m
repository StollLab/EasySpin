% trajectories of coupled stochastic rotational diffusion
%==========================================================================
% This script simulates two "levels" of rotational diffusion: fast,
% restricted local diffusion, and slow, Brownian global diffusion. The two
% models are coupled by performing quaternion multiplication between the
% sets of trajectories.

clear

% Local diffusion

Sys.tcorr = 1e-9;              % local rotational correlation time, s
Sys.Potential = [2,0,0,8.0];   % orientational potential with L=2, M=0, 
                               % K=0, and a coefficient of 2.0

Par.dt = Sys.tcorr/10;         % time step, s
Par.nSteps = 400;              % number of steps
Par.nTraj = 1;                 % number of trajectories
Par.OriStart = [0;0;0];        % starting orientation on the north pole

[t, RTrajLocal, qTrajLocal] = stochtraj_diffusion(Sys,Par);

% Global diffusion

Sys.tcorr = 20e-9;             % global rotational correlation time
Sys.Potential = [];            % remove the potential to simulate 
                               % unrestricted Brownian diffusion

[t, ~, qTrajGlobal] = stochtraj_diffusion(Sys,Par);


% combine trajectories using quaternion multiplication
qTraj = quatmult(qTrajGlobal,qTrajLocal);

% convert to trajectories of rotation matrices
RTraj = quat2rotmat(qTraj);

% extract the Z-axis vector of the body-fixed coordinate system
ZVecTrajLocal = squeeze(RTrajLocal(:,3,:,:));
ZVecTraj = squeeze(RTraj(:,3,:,:));

subplot(1,2,1)
plot3(ZVecTrajLocal(1,:),ZVecTrajLocal(2,:),ZVecTrajLocal(3,:));
axis equal
axlim = 1.2;
xlim([-1 1]*axlim);
ylim([-1 1]*axlim);
zlim([-1 1]*axlim);

subplot(1,2,2)
plot3(ZVecTraj(1,:),ZVecTraj(2,:),ZVecTraj(3,:));
axis equal
axlim = 1.2;
xlim([-1 1]*axlim);
ylim([-1 1]*axlim);
zlim([-1 1]*axlim);