% Propagate trajectories of orientational quaternions

% Implementation based on 
%  [1] Sezer, et al., J.Chem.Phys. 128, 165106 (2008)
%       http://dx.doi.org/10.1063/1.2908075
%  [2] Leimkuhler, Appl.Math.Res.Express 2013, 34 (2013)
%       https://doi.org/10.1093/amrx/abs010
%    newer method for SDE integration that starts with weak convergence of 
%    order 1, but approaches weak convergence of order 2 exponentially fast

function q = stochtraj_proprottraj(q, Sim, iter)

nSteps = Sim.nSteps;
nTraj = Sim.nTraj;
dt = Sim.dt;
Sim.Diff = Sim.Diff(:);
Diff = Sim.Diff;
lambda = Sim.lambda;
LMK = Sim.LMK;
if ~isfield(Sim,'interpGrad')
  interpGrad = [];
else
  interpGrad = Sim.interpGrad;
end

isEigenPot = ~isempty(lambda);
isNumericPot = ~isempty(interpGrad);

if iter>1
    % If propagation is being extended, initialize q from the last set
  startStep = 1;
else
  % First step has already been initialized by starting orientation, so
  % skip to second step
  startStep = 2;
end

%  Pre-calculate angular steps due to random torques
%  (Eq. 61 from reference, without factor of 1/2)
randAngSteps = genRandSteps(Sim);

for iStep = startStep:nSteps
  
  if iStep==1
    qLast = q(:,:,end);
    q = zeros(4,nTraj,nSteps);
  else
    qLast = q(:,:,iStep-1);
  end
  
  currRandAngStep = randAngSteps(:,:,iStep);
  
  % Calculate total torque
  if isEigenPot
    % use Wigner functions of quaternions to calculate torque
    torque = stochtraj_calcanistorque(LMK, lambda, qLast);
    AngStep = bsxfun(@times,torque,Diff*dt) + currRandAngStep;
  elseif isNumericPot
    % use orienting pseudopotential functions of Euler angles to calculate
    % torque
    [alpha, beta, gamma] = quat2euler(qLast,'active');
    pxint = interp3fast(interpGrad{1}, alpha, beta, gamma);
    pyint = interp3fast(interpGrad{2}, alpha, beta, gamma);
    pzint = interp3fast(interpGrad{3}, alpha, beta, gamma);
    torque = [-pxint.'; -pyint.'; -pzint.'];
    AngStep = bsxfun(@times,torque,Diff*dt) + currRandAngStep;
  else
    % If there is no orienting potential, then there is no torque to
    % calculate
    AngStep = currRandAngStep;
  end

  % Calculate size and normalized axis vector of angular step
  theta = sqrt(sum(AngStep.^2, 1));
  ux = AngStep(1,:,:)./theta;
  uy = AngStep(2,:,:)./theta;
  uz = AngStep(3,:,:)./theta;
  
  st = sin(theta/2);
  ct = cos(theta/2);

%   U = [    ct, -ux.*st, -uy.*st, -uz.*st; ...
%        ux.*st,      ct,  uz.*st, -uy.*st; ...
%        uy.*st, -uz.*st,      ct,  ux.*st; ...
%        uz.*st,  uy.*st, -ux.*st,      ct];

  % Calculate q for the first time step
  
  q1 = qLast(1,:);
  q2 = qLast(2,:);
  q3 = qLast(3,:);
  q4 = qLast(4,:);

  q(1,:,iStep) = q1.*ct - q2.*ux.*st - q3.*uy.*st - q4.*uz.*st;
  q(2,:,iStep) = q2.*ct + q1.*ux.*st - q4.*uy.*st + q3.*uz.*st;
  q(3,:,iStep) = q3.*ct + q4.*ux.*st + q1.*uy.*st - q2.*uz.*st;
  q(4,:,iStep) = q4.*ct - q3.*ux.*st + q2.*uy.*st + q1.*uz.*st;
  
%   diff = 1.0-sqrt(sum(q(:,:,iStep).*q(:,:,iStep), 1));
%   
%   if any(abs(diff(:)) > 1e-14)
%     error('Quaternions are not normalized on step %d.\n', iStep)
%   end
  
end

end

function randAngStep = genRandSteps(Sim)
% generate random angular steps using one of two different MC integrators

% generate Gaussian random deviates
randns = randn(3,Sim.nTraj,Sim.nSteps);

switch Sim.Integrator
  case 'Euler-Maruyama'
    % standard method for integrating a first order SDE
    randAngStep = bsxfun(@times, randns, sqrt(2*Sim.Diff*Sim.dt));
  case 'Leimkuhler-Matthews'
    % see Ref. [2]
    randAngStep = bsxfun(@times, (randns(:,:,1:end-1)+randns(:,:,2:end))/2, ...
                                     sqrt(2*Sim.Diff*Sim.dt));
end

end

function Vq = interp3fast(F, Xq, Yq, Zq)
% Adapted from MATLAB's interp3 function, and makes the following
% assumptions regarding input:
%    Monotonic grid vectors X, Y, Z were fed to griddedInterpolant to
%    obtain F
%    V is ndgrid-ordered, not meshgrid-ordered (fed to griddedInterpolant)
%interp3fast 3-D interpolation (table lookup).

Vq = F([Xq.',Yq.',Zq.']);

end

