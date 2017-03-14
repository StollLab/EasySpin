% propagate_quantum  Propagate the density matrix of a spin-1/2 14N
%                    nitroxide using different methods from the literature.
%
%   rho_t = propagate_quantum(Sys,Par,Opt,omega,RTraj);
%
%     omega          double
%                    microwave frequency for CW field sweep, in Hz
%
%     CenterField    double
%                    magnetic field, in mT
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
%                    'Sezer': propagate using the m_s=-1/2 subspace
%                    'DeSensi': propagate using an eigenvalue method
%                    'Oganesyan': propagate using correlation functions
%
%
%
%   Output:
%     rho_t          numeric, size = (3,3,nTraj,nSteps)
%                    a series of density matrices

% Implementations based on 
% [1] Sezer, et al., J. Chem. Phys. 128, 165106 (2008)
%      http://dx.doi.org/10.1063/1.2908075 
% [2] DeSensi, et al., Biophys. J. 94, 3798 (2008)
%      http://dx.doi.org/10.1529/biophysj.107.125419
% [3] Oganesyan, Phys. Chem. Chem. Phys. 13, 4724 (2011)
%      http://dx.doi.org/10.1039/c0cp01068e

function rho_t = propagate_quantum(Sys, Par, Opt, omega, CenterField, RTraj)

% Preprocessing
% -------------------------------------------------------------------------

if nargin~=6
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

g = Sys.g;
A = Sys.A;

if ~isfield(Par,'dt'), error('Time step not specified.'); end

dt = Par.dt;
nTraj = Par.nTraj;
nSteps = Par.nSteps;

if ~isequal(size(g),[1,3]) || ~isequal(size(A),[1,3])
  error('g and A tensor values must be 3-vectors.')
end

RTrajInv = permute(RTraj,[2,1,3,4]);

% Gamma = gfree*bmagn/(planck/2/pi);  % rad s^-1 T^-1
% omegaN = 19.331e6*B;  % gyromagnetic ratio for 14N: 
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
                                         RTrajInv))*1e6*2*pi;  % MHz (s^-1) -> Hz (rad s^-1)

    Gp_zz = Gp_tensor(3,3,:,:);

    a = sqrt(A_tensor(1,3,:,:).*A_tensor(1,3,:,:) ...
           + A_tensor(2,3,:,:).*A_tensor(2,3,:,:) ...
           + A_tensor(3,3,:,:).*A_tensor(3,3,:,:));

%     theta = Gamma*dt*0.5*squeeze(a);  % use this if A is given in G
    theta = dt*0.5*squeeze(a);
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


    % Calculate propagator
%     U = bsxfun(@times, exp(-1i*dt*0.5*Gamma*B*Gp_zz), exp_array);
      U = bsxfun(@times, exp(-1i*dt*0.5*omega*Gp_zz), exp_array);
%       U2 = matmult(U,U);

    % Eq. 34 in [1]
    
    for iStep=2:nSteps
     rho_t(:,:,:,iStep) = mmult(U(:,:,:,iStep-1),...
                                mmult(rho_t(:,:,:,iStep-1), U(:,:,:,iStep-1)));

    %  Trace of a product of matrices is the sum of entry-wise products
%      rho_t(:,:,:,iStep) = U2(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
    %  Full sandwich product
    %  rho_t(:,:,:,iStep) = matmult(U(:,:,:,iStep-1), ...
    %                                   matmult(rho_t(:,:,:,iStep-1), ...
    %                                           U(:,:,:,iStep-1)));

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'DeSensi'  % see Ref [2]

    rho_t = zeros(6,6,nTraj,nSteps);
    rho_t(1:3,4:6,:,1) = repmat(eye(3),1,1,nTraj,1);
    rho_t(4:6,1:3,:,1) = repmat(eye(3),1,1,nTraj,1);
%     rho_t = zeros(3,3,nTraj,nSteps);
%     rho_t(:,:,:,1) = repmat(eye(3),1,1,nTraj,1);

    g_tr = sum(g);

    % Perform rotations on g- and A-tensors
%     g_tensor = matmult(RTraj, matmult(repmat(diag(g),1,1,nTraj,nSteps), ...
%                                       RTrajInv));
    Gp_tensor = matmult(RTraj, matmult(repmat(diag(g),1,1,nTraj,nSteps), ...
                                      RTrajInv))/gfree - g_tr/3/gfree;

    A_tensor = matmult(RTraj, matmult(repmat(diag(A),1,1,nTraj,nSteps), ...
                                      RTrajInv))*1e6*2*pi;  % MHz (s^-1) -> Hz (rad s^-1)

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

    ca = A_zz;% - 4*ma*omegaN;
    cb = A_zz;% - 4*mb*omegaN;

    ella = sqrt( (A_xz).^2 + (A_yz).^2 + ca.^2 );
    ellb = sqrt( (A_xz).^2 + (A_yz).^2 + cb.^2 );

    % Eigenvalues
%     Lambda  = [ 0.5*(geff + ellp);          0;                 0;                  0;           0;                  0;
%                                 0; 0.5*(geff);                 0;                  0;           0;                  0;
%                                 0;          0; 0.5*(geff + ellm);                  0;           0;                  0;
%                                 0;          0;                 0; -0.5*(geff + ellp);           0;                  0;
%                                 0;          0;                 0;                  0; -0.5*(geff);                  0;
%                                 0;          0;                 0;                  0;           0; -0.5*(geff + ellm) ];
    Lambda = zeros(6,6,nTraj,nSteps);
    Lambda(1,1,:,:) = exp(1i*dt*ma*(omega*geff + mp*ella));
    Lambda(2,2,:,:) = exp(1i*dt*ma*(omega*geff + m0*ella));
    Lambda(3,3,:,:) = exp(1i*dt*ma*(omega*geff + mm*ella));
%     Lambda = zeros(3,3,nTraj,nSteps);
    Lambda(4,4,:,:) = exp(1i*dt*mb*(omega*geff + mp*ellb));
    Lambda(5,5,:,:) = exp(1i*dt*mb*(omega*geff + m0*ellb));
    Lambda(6,6,:,:) = exp(1i*dt*mb*(omega*geff + mm*ellb));

    b   = (A_xz + 1i*A_yz);
    bst = (A_xz - 1i*A_yz);

    % Matrix of eigenvectors
    V  = [ (ca+ella).^2./b.^2,                     -bst./b,      (ca-ella).^2./b.^2, zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps);
           sqrt(2)*(ca+ella)./b,             sqrt(2)*ca./b,    sqrt(2)*(ca-ella)./b, zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps);
           ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps);
          zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps),      (cb+ellb).^2./b.^2,                 -bst./b,      (cb-ellb).^2./b.^2;
          zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps),    sqrt(2)*(cb+ellb)./b,           sqrt(2)*cb./b,    sqrt(2)*(cb-ellb)./b;
          zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),   ones(1,1,nTraj,nSteps) ];
%     V  = [    (cb+ellb).^2./b.^2,                 -bst./b,      (cb-ellb).^2./b.^2;
%             sqrt(2)*(cb+ellb)./b,           sqrt(2)*cb./b,    sqrt(2)*(cb-ellb)./b;
%           ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),   ones(1,1,nTraj,nSteps) ];

    % Normalize eigenvectors
    V = bsxfun(@rdivide,V,sqrt(sum(V.*conj(V),1)));
    
    % Calculate propagator
    Q = mmult(V, mmult(Lambda, conj(permute(V,[2,1,3,4]))));
    
    % Eq. 34 in reference
% FIXME round-off error might propagate through the matrix multiplication, need some sort of error control
    for iStep=2:nSteps
%       for iTraj=1:nTraj
      rho_t(:,:,:,iStep) = mmult(Q(:,:,:,iStep-1),...
                                   mmult(rho_t(:,:,:,iStep-1),...
                                         conj(permute(Q(:,:,:,iStep-1),[2,1,3,4]))));
%         Q = V(:,:,iTraj,iStep-1)*Lambda(:,:,iTraj,iStep-1)*V(:,:,iTraj,iStep-1)';
%         rho_t(:,:,iTraj,iStep) = Q*rho_t(:,:,iTraj,iStep-1)*Q;                    
%       end
    end
    
    rho_t = rho_t(4:6,1:3,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'Oganesyan'  % see Ref [3]
    rho_t = zeros(6,6,nTraj,nSteps);
    rho_t(1:3,4:6,:,1) = repmat(eye(3),1,1,nTraj,1);
    rho_t(4:6,1:3,:,1) = repmat(eye(3),1,1,nTraj,1);
    rho_t = rho_t(:);
    
    B0 = {0,0,Centerfield/1e3};  % mT -> T
    
    SpinOps{1,1} = sop([1/2,1],'xe');
    SpinOps{1,2} = sop([1/2,1],'ye');
    SpinOps{1,3} = sop([1/2,1],'ze');
    
    IncludeNuclearZeeman = 0;
    ExplicitFieldSweep = 0;
    
    [T0,T1,T2,F0,F1,F2,isFieldDep] = magint(Sys,SpinOps,CenterField, ...
                                            IncludeNuclearZeeman, ...
                                            ExplicitFieldSweep);
    
    

end