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
%   Opt: optional settings
%     Method         string
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

function rho_t = propagate_quantum(Sys, Par, Opt, RTraj)

% Preprocessing
% -------------------------------------------------------------------------

if nargin~=4
  error('Wrong number of input arguments.')
end

Method = Opt.Method;

% Check shapes of inputs
if size(RTraj,1)~=3 || size(RTraj,2)~=3 || size(RTraj,3)~=Par.nTraj ...
    || size(RTraj,4)~=Par.nSteps
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
% Gamma = 1.7608597e11;
omega0 = mt2mhz(B*1e3)*1e6;
omegaN = 19.331e6*B;  % gyromagnetic ratio for 14N: 
                      % 19.331x10^6 rad s^-1 T^-1


% Simulation
% -------------------------------------------------------------------------

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

    % Note: here and below, we use 2*dt instead of dt for the time step 
    % since M_+ = tr(U*rho*U) = tr(U*U*rho) = tr(U^2*rho)
    % This is follows from the fact that the trace of a product of matrices 
    % is invariant under cyclic permutations
    theta = Gamma*dt*0.5*squeeze(a);
    nx = squeeze(A_tensor(1,3,:,:)./a);
    ny = squeeze(A_tensor(2,3,:,:)./a);
    nz = squeeze(A_tensor(3,3,:,:)./a);

    % Eqs. A1-A2 in [1]
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
    U = bsxfun(@times, exp(-1i*dt*0.5*Gamma*B*Gp_zz), exp_array);


    % Eq. 34 in [1]
    
    for iStep=2:nSteps
    %  Use of permutation property of a trace of a product
     rho_t(:,:,:,iStep) = matmult(U(:,:,:,iStep-1), matmult(rho_t(:,:,:,iStep-1), U(:,:,:,iStep-1)));

    %  Trace of a product of matrices is the sum of entry-wise products
    %  rho_t(:,:,:,iStep) = U(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
    %  Full sandwich product
    %  rho_t(:,:,:,iStep) = matmult(U(:,:,:,iStep-1), ...
    %                                   matmult(rho_t(:,:,:,iStep-1), ...
    %                                           U(:,:,:,iStep-1)));

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'DeSensi'  % see Ref [2]

%     rho_t = zeros(6,6,nTraj,nSteps);
%     rho_t(1:3,4:6,:,1) = repmat(eye(3),1,1,nTraj,1);
%     rho_t(4:6,1:3,:,1) = repmat(eye(3),1,1,nTraj,1);
    rho_t = zeros(3,3,nTraj,nSteps);
    rho_t(:,:,:,1) = repmat(eye(3),1,1,nTraj,1);

    g_tr = sum(g);

    % Perform rotations on g- and A-tensors
%     g_tensor = matmult(RTraj, matmult(repmat(diag(g),1,1,nTraj,nSteps), ...
%                                       RTrajInv));
    Gp_tensor = matmult(RTraj, matmult(repmat(diag(g),1,1,nTraj,nSteps), ...
                                      RTrajInv))/gfree - g_tr/3/gfree;

    A_tensor = matmult(RTraj, matmult(repmat(diag(A),1,1,nTraj,nSteps), ...
                                      RTrajInv));

    Gp_zz = Gp_tensor(3,3,:,:);
%     g_zz = g_tensor(3,3,:,:);

    A_xz = A_tensor(1,3,:,:);
    A_yz = A_tensor(2,3,:,:);
    A_zz = A_tensor(3,3,:,:);
    
%     SzIe = sop([1/2,1], 'ze');
%     SzIx = sop([1/2,1], 'zx');
%     SzIy = sop([1/2,1], 'zy');
%     SzIz = sop([1/2,1], 'zz');
%     SeIz = sop([1/2,1], 'ez');

    % Eq. 22 in [2]
%     H = bsxfun(@times,(g_zz*bmagn/planck/2/pi*B-omega0), SzIe) ...
%     H = bsxfun(@times,Gamma*Gp_zz*B, SzIe) ...
%         + bsxfun(@times,Gamma*A_xz, SzIx) ...
%         + bsxfun(@times,Gamma*A_yz, SzIy) ...
%         + bsxfun(@times,Gamma*A_zz, SzIz);
%         - bsxfun(@times,omegaN, SeIz);

    % Eq. 24-27 in [2]
%     geff = -bmagn/(planck/2/pi)*B*g_zz + omega0; 
    geff = Gp_zz;
    
    ma = 0.5;
    mb = -0.5;
    
    mp = 1;
    m0 = 0;
    mm = -1;

    ca = Gamma*A_zz; % - 4*ma*omegaN;
    cb = Gamma*A_zz; % - 4*mb*omegaN;

    ella = sqrt( (Gamma*A_xz).^2 + (Gamma*A_yz).^2 + ca.^2 );
    ellb = sqrt( (Gamma*A_xz).^2 + (Gamma*A_yz).^2 + cb.^2 );

    % Eigenvalues
%     Lambda  = [ 0.5*(geff + ellp);          0;                 0;                  0;           0;                  0;
%                                 0; 0.5*(geff);                 0;                  0;           0;                  0;
%                                 0;          0; 0.5*(geff + ellm);                  0;           0;                  0;
%                                 0;          0;                 0; -0.5*(geff + ellp);           0;                  0;
%                                 0;          0;                 0;                  0; -0.5*(geff);                  0;
%                                 0;          0;                 0;                  0;           0; -0.5*(geff + ellm) ];
%     Lambda = zeros(6,6,nTraj,nSteps);
    Lambda = zeros(3,3,nTraj,nSteps);
%     dt=2*dt;
%     Lambda(1,1,:,:) = exp(-1i*dt*ma*(geff + mp*ella));
%     Lambda(2,2,:,:) = exp(-1i*dt*ma*(geff + m0*ella));
%     Lambda(3,3,:,:) = exp(-1i*dt*ma*(geff + mm*ella));
    Lambda(1,1,:,:) = exp(1i*dt*mb*(Gamma*B*geff + mp*ellb));
    Lambda(2,2,:,:) = exp(1i*dt*mb*(Gamma*B*geff + m0*ellb));
    Lambda(3,3,:,:) = exp(1i*dt*mb*(Gamma*B*geff + mm*ellb));

    b   = Gamma*(A_xz + 1i*A_yz);
    bst = Gamma*(A_xz - 1i*A_yz);

    % Matrix of eigenvectors
%     V  = [ (ca+ella).^2./b.^2,                     -bst./b,      (ca-ella).^2./b.^2, zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps);
%            sqrt(2)*(ca+ella)./b,             sqrt(2)*ca./b,    sqrt(2)*(ca-ella)./b, zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps);
%            ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps);
%           zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps),      (cb+ellb).^2./b.^2,                 -bst./b,      (cb-ellb).^2./b.^2;
%           zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps),    sqrt(2)*(cb+ellb)./b,           sqrt(2)*cb./b,    sqrt(2)*(cb-ellb)./b;
%           zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),   ones(1,1,nTraj,nSteps) ];
    V  = [    (cb+ellb).^2./b.^2,                 -bst./b,      (cb-ellb).^2./b.^2;
            sqrt(2)*(cb+ellb)./b,           sqrt(2)*cb./b,    sqrt(2)*(cb-ellb)./b;
          ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),   ones(1,1,nTraj,nSteps) ];
 


%     % Matrix of eigenvector norms
%     norm = zeros(6,6,nTraj,nSteps);
%     norm = [ abs((ca+ella).^2./b.^2 + sqrt(2)*(ca+ella)./b + ones(1,1,nTraj,nSteps)),...
%                                abs(-bst./b + sqrt(2)*ca./b + ones(1,1,nTraj,nSteps)),...
%              abs((ca-ella).^2./b.^2 + sqrt(2)*(ca-ella)./b + ones(1,1,nTraj,nSteps)),...
%              abs((cb+ellb).^2./b.^2 + sqrt(2)*(cb+ellb)./b + ones(1,1,nTraj,nSteps)),...
%                                abs(-bst./b + sqrt(2)*cb./b + ones(1,1,nTraj,nSteps)),...
%              abs((cb-ellb).^2./b.^2 + sqrt(2)*(cb-ellb)./b + ones(1,1,nTraj,nSteps)) ];
        
    V = bsxfun(@rdivide,V,sqrt(sum(V.*conj(V),1)));
%     V = bsxfun(@rdivide,V,norm);
    
%     [v,D] = eig(H(:,:,1,1));
%     
%     (V(:,:,1,1)-v)./V(:,:,1,1)
%     (Lambda(:,:,1,1)-D)./Lambda(:,:,1,1)
    
%     V_test = zeros(6,6);
%     
%     for iStep=1:nSteps
%       for iTraj=1:nTraj
%         V_test(:,:,iTraj,iStep) = V(:,:,iTraj,iStep)*V(:,:,iTraj,iStep)' - eye(6);
%       end
%     end
% 
%     if any(V_test(:) > 1e-10)
%       error('V is not unitary.')
%     end

    % Calculate propagator
%     Q = zeros(6,6,nTraj,nSteps);
%     for iStep=1:nSteps
%       for iTraj=1:nTraj
%         Q(:,:,iTraj,iStep) = V(:,:,iTraj,iStep)*Lambda(:,:,iTraj,iStep)*V(:,:,iTraj,iStep)';
%       end
%     end
    
%     Q_test = zeros(6,6,nTraj,nSteps);
    
%     for iStep=1:nSteps
%       for iTraj=1:nTraj
%         Q_test(:,:,iTraj,iStep) = Q(:,:,iTraj,iStep)*Q(:,:,iTraj,iStep)' - eye(6);
%       end
%     end
                 
%     if any(Q_test(:) > 1e-10)
%       error('Q is not unitary.')
%     end

%     V = V(4:6,4:6,:,:);
%     Lambda = Lambda(4:6,4:6,:,:);
    
    Q = matmult(V, matmult(Lambda, conj(permute(V,[2,1,3,4]))));
    % Eq. 34 in reference
% FIXME round-off error propagates through the matrix multiplication, need some sort of error control
    for iStep=2:nSteps
%       for iTraj=1:nTraj
%         Q = expm(-1i*dt*H(4:6,4:6,iTraj,iStep-1));
        % Note that that we need Q*rho*Q, not Q*rho*Q', if are working in the decoupled m_beta subspace
%         rho_t(:,:,:,iStep) = matmult(Q(:,:,:,iStep),rho_t(:,:,:,iStep-1));
      rho_t(:,:,:,iStep) = matmult(Q(:,:,:,iStep-1),matmult(rho_t(:,:,:,iStep-1),Q(:,:,:,iStep-1)));
%         rho_t(:,:,iTraj,iStep) = Q*rho_t(:,:,iTraj,iStep-1)*Q;
%         Q = V(:,:,iTraj,iStep)*Lambda(:,:,iTraj,iStep)*V(:,:,iTraj,iStep)';
%         rho_t(:,:,iTraj,iStep) = Q*rho_t(:,:,iTraj,iStep-1)*Q';
%       end
    end
    
%     rho_t = rho_t(4:6,1:3,:,:);
    
end

end