% propagate_quantum  Propagate the density matrix of a spin-1/2 14N
%                    nitroxide using different methods from the literature.
%
%   rho_t = propagate_quantum(Sys,Par,RTraj,Method);
%
%     RTraj          numeric, size = (3,3,nTraj,nSteps)
%                    a series of rotation matrices
%
%   Sys: stucture with system's dynamical parameters
%     g              numeric, size = (1,3)
%                    principal values of the g-tensor
%     A              numeric, size = (1,3)
%                    principal values of the A-tensor
%
%   Par: structure with simulation parameters
%     dt             double
%                    time step (in seconds)
%
%   Exp: experimental parameter settings
%     B              double  TODO can this be replaced by a fieldsweep and then used to extract omega0?
%                    center magnetic field
%
%   Method           string
%                    'Sezer': propagate the density matrix using the 
%                     m_s=-1/2 subspace
%                    'DeSensi': propagate the density matrix using an 
%                     eigenvalue method
%
%
%
%   Output:
%     rho_t          numeric, size = (3,3,nTraj,nSteps)
%                    a series of density matrices

% Implementation based on 
% [1] Sezer, et al., J.Chem.Phys. 128, 165106 (2008)
%      http://dx.doi.org/10.1063/1.2908075 
% [2] DeSensi, et al., Biophys. J. 94, 3798 (2008)
%      http://dx.doi.org/10.1529/biophysj.107.125419

function rho_t = propagate_quantum(Sys, Par, RTraj, Method)
%% Preprocessing
%========================================================================

if nargin~=4
  error('Wrong number of input arguments.')
end

% Check shapes of inputs
if size(RTraj,1)~=3 || size(RTraj,2)~=3 || size(RTraj,3)~=Par.nTraj || size(RTraj,4)~=Par.nSteps
  error('Array of rotation matrices must be of size = (3,3,nTraj,nSteps).')
end

if ~isfield(Sys,'g'), error('g-tensor not specified.'); end
if ~isfield(Sys,'A'), error('A-tensor not specified.'); end
if ~isfield(Sys,'B'), error('Magnetic field not specified.'); end

g = Sys.g;
A = Sys.A;
B = Sys.B;

if ~isfield(Par,'dt'), error('Time step not specified.'); end

dt = Par.dt;
nTraj = Par.nTraj;
nSteps = Par.nSteps;

if ~isequal(size(g),[1,3]) || ~isequal(size(A),[1,3])
  error('g and A tensor values must be 3-vectors.')
end

% Check for orthogonality of rotation matrices
RTrajInv = permute(RTraj,[2,1,3,4]);

rot_mat_test = matmult(RTraj,RTrajInv) ...
               - repmat(eye(3),1,1,nTraj,nSteps);

if any(rot_mat_test > 1e-10)
  error('The rotation matrices are not orthogonal.')
end


Gamma = gfree*bmagn/(planck/2/pi);  % rad s^-1 T^-1

%% Simulation
%========================================================================

switch Method
  case 'Sezer'  % see Ref [1]

    rho_t = zeros(3,3,nTraj,nSteps);
    rho_t(:,:,:,1) = repmat(eye(3),1,1,nTraj,1);

    g_tr = sum(g);

    % Perform rotations on g- and A-tensors
    Gp_tensor = matmult(RTraj, matmult(repmat(diag(g),1,1,nTraj,nSteps), ...
                                          RTrajInv))/gfree - g_tr/3/gfree;

    A_tensor = matmult(RTraj, matmult(repmat(diag(A),1,1,nTraj,nSteps), ...
                                         RTrajInv));

    Gp_zz = Gp_tensor(3,3,:,:);

    a = sqrt(A_tensor(1,3,:,:).*A_tensor(1,3,:,:) ...
           + A_tensor(2,3,:,:).*A_tensor(2,3,:,:) ...
           + A_tensor(3,3,:,:).*A_tensor(3,3,:,:));

    % Note: here and below, we use 2*dt instead of dt for the time step since
    % M_+ = tr(U*rho*U) = tr(U*U*rho) = tr(U^2*rho)
    % This is follows from the fact that the trace of a product of matrices is 
    % invariant under cyclic permutations
    theta = Gamma*(2*dt)*0.5*squeeze(a);
    nx = squeeze(A_tensor(1,3,:,:)./a);
    ny = squeeze(A_tensor(2,3,:,:)./a);
    nz = squeeze(A_tensor(3,3,:,:)./a);

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


    % Again, here we use 2*dt instead of dt
    U = bsxfun(@times, exp(-1i*(2*dt)*0.5*Gamma*B*Gp_zz), exp_array);


    % Eq. 34 in reference
    for iStep=2:nSteps
    %  Use of permutation property of a trace of a product
     rho_t(:,:,:,iStep) = matmult(U(:,:,:,iStep-1),rho_t(:,:,:,iStep-1));

    %  Trace of a product of matrices is the sum of entry-wise products
    %  rho_t(:,:,:,iStep) = U(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
    %  Full sandwich product
    %  rho_t(:,:,:,iStep) = matmult(U(:,:,:,iStep-1), ...
    %                                   matmult(rho_t(:,:,:,iStep-1), ...
    %                                           U(:,:,:,iStep-1)));

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'DeSensi'  % see Ref [2]
    error('DeSensi method not implemented yet!')

    rho_t = zeros(6,6,nTraj,nSteps);
    rho_t(1:3,4:6,:,1) = repmat(eye(3),1,1,nTraj,1);
    rho_t(4:6,1:3,:,1) = repmat(eye(3),1,1,nTraj,1);

    % g_tr = sum(g);

    % Perform rotations on g- and A-tensors
    g_tensor = matmult(RTraj, matmult(repmat(diag(g),1,1,nTraj,nSteps), ...
                                      RTrajInv));

    A_tensor = matmult(RTraj, matmult(repmat(diag(A),1,1,nTraj,nSteps), ...
                                      RTrajInv));

    g_zz = g_tensor(3,3,:,:);

    A_xz = A_tensor(1,3,:,:);
    A_yz = A_tensor(2,3,:,:);
    A_zz = A_tensor(3,3,:,:);

    % Eq. 22 in [2]
    H = (g_zz*bmagn*B - omega0)*S_z + Gamma*A_xz*I_x*S_z + Gamma*A_yz*I_y*S_z ...
                                   + Gamma*A_zz*I_z*S_z - omegaN*I_z;

    % Eq. 24-27 in [2]
    geff = bmagn*B*g_zz - omega0;

    cp = Gamma*A_zz - 4*(0.5)*omegaN;
    cm = Gamma*A_zz - 4*(-0.5)*omegaN;

    ellp = sqrt( (Gamma*A_xz).^2 + (Gamma*A_yz).^2 + cp.^2 );
    ellm = sqrt( (Gamma*A_xz).^2 + (Gamma*A_yz).^2 + cm.^2 );

    % Eigenvalues
    lambda(1,1,:,:) =  0.5*(geff + ellp);
    lambda(2,1,:,:) =  0.5*(geff);
    lambda(3,1,:,:) =  0.5*(geff + ellm);
    lambda(4,1,:,:) = -0.5*(geff + ellp);
    lambda(5,1,:,:) = -0.5*(geff);
    lambda(6,1,:,:) = -0.5*(geff + ellm);

    b   = A_xz + 1i*A_yz;
    bst = A_xz - 1i*A_yz;

    % Eigenvectors
    nu = [ (cp+lp).^2./(4*b.^2),         -bst./b, (cp-lp).^2./(4*b.^2),                    0,               0,                   0;
           (cp+lp)./(sqrt(2)*b), cp./(sqrt(2)*b), (cp-lp)./(sqrt(2)*b),                    0,               0,                   0;
                              1,               1,                    1,                    0,               0,                   0;
                              0,               0,                    0, (cm+lm).^2./(4*b.^2),         -bst./b, (cm-lm).^2./(4*b.^2);
                              0,               0,                    0, (cm+lm)./(sqrt(2)*b), cm./(sqrt(2)*b), (cm-lm)./(sqrt(2)*b);
                              0,               0,                    0,                    1,               1,                   1 ];


    % Again, here we use 2*dt instead of dt
    U = bsxfun(@times, exp(-1i*(2*dt)*0.5*Gamma*B*Gp_zz), exp_array);


    % Eq. 34 in reference
    for iStep=2:nSteps
    %  Use of permutation property of a trace of a product
     rho_t(:,:,:,iStep) = matmult(U(:,:,:,iStep-1),rho_t(:,:,:,iStep-1));

    %  Trace of a product of matrices is the sum of entry-wise products
    %  rho_t(:,:,:,iStep) = U(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
    %  Full sandwich product
    %  rho_t(:,:,:,iStep) = matmult(U(:,:,:,iStep-1), ...
    %                                   matmult(rho_t(:,:,:,iStep-1), ...
    %                                           U(:,:,:,iStep-1)));

    end

end

end