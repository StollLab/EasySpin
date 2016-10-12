% propagate_hfield  Propagate the density matrix for a spin-1/2 nitroxide
%                   using the high field approximation.
%
%  rho_t = propagate_hfield(Sys, RTraj);
%
%  Sys: structure with simulation parameters
%     g              1x3 array, principal values of the g-tensor
%     A              1x3 array, principal values of the A-tensor
%     B              douvle, magnetic field
%     dt             double, time step (in seconds)
%     RTraj          3x3xnStepsxnSims array, a series of rotation matrices
%
%  Output:
%     rho_t          a series of density matrices
%
% Implementation based on 
% - Sezer, Freed, Roux, J.Chem.Phys. 128, 165106 (2008)

function rho_t = propagate_hfield(Sys, RTraj)

% Check shapes of inputs
if (size(RTraj,1)~=size(RTraj,2)) || size(RTraj,1)~=3
  error('Size of rotation matrices must be 3x3xnStepsxnTraj.')
end

nSteps = size(RTraj,3);
nTraj = size(RTraj,4);

if ~isfield(Sys,'g'), error('g-tensor not specified.'); end
if ~isfield(Sys,'A'), error('A-tensor not specified.'); end
if ~isfield(Sys,'B'), error('Magnetic field not specified.'); end
if ~isfield(Sys,'dt'), error('Time step not specified.'); end

g = Sys.g;
A = Sys.A;

if ~isequal(size(g),[1,3]) || ~isequal(size(A),[1,3])
    error('g and A tensor values must be 3-vectors.')
end

% Check for orthogonality of rotation matrices
RTrajInv = permute(RTraj,[2,1,3,4]);

rot_mat_test = matmult(RTraj,RTrajInv) ...
               - repmat(eye(3),1,1,nSteps,nTraj);

if any(rot_mat_test > 1e-10);
    error('The rotation matrices are not orthogonal.')
end

dt = Sys.dt;
B = Sys.B;

%========================================================================
% Begin calculation
%========================================================================

Gamma = 1.7608597e11; % rad s^-1 T^-1
ge = 2.0023193;

rho_t = zeros(3,3,nSteps,nTraj);
rho_t(:,:,1,:) = repmat(eye(3),1,1,1,nTraj);

g_tr = sum(g);

% Perform rotations on g- and A-tensors
Gp_tensor = matmult(RTraj, matmult(repmat(diag(g),1,1,nSteps,nTraj), ...
                                      RTrajInv))/ge - g_tr/3/ge;

A_tensor = matmult(RTraj, matmult(repmat(diag(A),1,1,nSteps,nTraj), ...
                                     RTrajInv));

Gp_zz = Gp_tensor(3,3,:,:);

a = sqrt(A_tensor(1,3,:,:).*A_tensor(1,3,:,:) ...
       + A_tensor(2,3,:,:).*A_tensor(2,3,:,:) ...
       + A_tensor(3,3,:,:).*A_tensor(3,3,:,:));

theta = Gamma*dt*0.5*squeeze(a);
nx = squeeze(A_tensor(1,3,:,:)./a);
ny = squeeze(A_tensor(2,3,:,:)./a);
nz = squeeze(A_tensor(3,3,:,:)./a);

% Use 2*theta to obtain U^2 = exp[-i*(2*theta)*N]
% Eqs. A1-A2 in reference
ct = cos(theta) - 1;
st = -sin(theta);

exp_array(1,1,:,:) = 1 + ct.*(nz.*nz + 0.5*(nx.*nx + ny.*ny)) ...
                   + 1i*st.*nz;
exp_array(1,2,:,:) = sqrt(0.5)*(st.*ny + ct.*nz.*nx) ...
                   + 1i*sqrt(0.5)*(st.*nx - ct.*nz.*ny);
exp_array(1,3,:,:) = 0.5*ct.*(nx.*nx - ny.*ny) ... 
                   - 1i*ct.*nx.*ny;
exp_array(2,1,:,:) = sqrt(0.5)*(-st.*ny + ct.*nz.*nx) ...
                   + 1i*sqrt(0.5)*(st.*nx + ct.*nz.*ny);
exp_array(2,2,:,:) = 1 + ct.*(nx.*nx + ny.*ny);
exp_array(2,3,:,:) = sqrt(0.5)*(st.*ny - ct.*nz.*nx) ...
                   + 1i*sqrt(0.5)*(st.*nx + ct.*nz.*ny);
exp_array(3,1,:,:) = 0.5*ct.*(nx.*nx - ny.*ny) ...
                   + 1i*ct.*nx.*ny;
exp_array(3,2,:,:) = sqrt(0.5)*(-st.*ny - ct.*nz.*nx) ...
                   + 1i*sqrt(0.5)*(st.*nx - ct.*nz.*ny);
exp_array(3,3,:,:) = 1 + ct.*(nz.*nz + 0.5*(nx.*nx + ny.*ny))  ...
                   - 1i*st.*nz;


% Since tr(U*rho*U) = tr(U^2*rho), U^2 is calcuated by using
% 2*theta above and 2*dt below
% U2 = bsxfun(@times, exp(-1i*2*dt*0.5*Gamma*B*Gp_zz), exp_array);

U = bsxfun(@times, exp(-1i*dt*0.5*Gamma*B*Gp_zz), exp_array);


% Eq. 34 in reference
for iStep=2:nSteps
%   rho_t(:,:,iStep,:) = U2(:,:,iStep-1,:).*rho_t(:,:,iStep-1,:);
%   rho_t(:,:,iStep,:) = matmult(U2(:,:,iStep-1,:),rho_t(:,:,iStep-1,:));
  rho_t(:,:,iStep,:) = matmult(U(:,:,iStep-1,:), ...
                                  matmult(rho_t(:,:,iStep-1,:), ...
                                          U(:,:,iStep-1,:)));

end

end