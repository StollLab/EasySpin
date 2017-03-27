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
%     qTraj          numeric, size = (4,nTraj,nSteps)
%                    a series of quaternions representing orientations
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

function rho_t = propagate_quantum(Sys, Par, Opt, omega, CenterField)

% Preprocessing
% -------------------------------------------------------------------------

if nargin~=5
  error('Wrong number of input arguments.')
end

Method = Opt.Method;

if isfield(Par,'RTraj')
  RTraj = Par.RTraj;
  RTrajInv = permute(RTraj,[2,1,3,4]);
  if size(RTraj,1)~=3 || size(RTraj,2)~=3 || size(RTraj,3)~=Par.nTraj ...
      || size(RTraj,4)~=Par.nSteps
    error('Array of rotation matrices must be of size = (3,3,nTraj,nSteps).')
  end
end
if isfield(Par,'qTraj')
  qTraj = Par.qTraj;
  
  if size(qTraj,1)~=4 || size(qTraj,2)~=Par.nTraj ...
      || size(qTraj,3)~=Par.nSteps
    error('Array of quaternions must be of size = (4,nTraj,nSteps).')
  end
else
  error('Either Par.RTraj or Par.qTraj must be provided.')
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
                                 mmult(rho_t(:,:,:,iStep-1), U(:,:,:,iStep-1),'complex'),...
                                 'complex');

    %  Trace of a product of matrices is the sum of entry-wise products
%      rho_t(:,:,:,iStep) = U2(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'DeSensi'  % see Ref [2]

    rho_t = zeros(6,6,nTraj,nSteps);
    rho_t(1:3,4:6,:,1) = repmat(eye(3),1,1,nTraj,1);
    rho_t(4:6,1:3,:,1) = repmat(eye(3),1,1,nTraj,1);

    g_tr = sum(g);

    % Perform rotations on g- and A-tensors
    Gp_tensor = matmult(RTraj, matmult(repmat(diag(g),1,1,nTraj,nSteps), ...
                                      RTrajInv))/gfree - g_tr/3/gfree;

    A_tensor = matmult(RTraj, matmult(repmat(diag(A),1,1,nTraj,nSteps), ...
                                      RTrajInv))*1e6*2*pi;  % MHz (s^-1) -> rad s^-1

    Gp_zz = Gp_tensor(3,3,:,:);
%     g_zz = g_tensor(3,3,:,:);

    A_xz = A_tensor(1,3,:,:);
    A_yz = A_tensor(2,3,:,:);
    A_zz = A_tensor(3,3,:,:);

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

    % Normalize eigenvectors
    V = bsxfun(@rdivide,V,sqrt(sum(V.*conj(V),1)));
    
    % Calculate propagator
    Q = mmult(V, mmult(Lambda, conj(permute(V,[2,1,3,4])), 'complex'), 'complex');
    
    % Eq. 34 in reference
% FIXME round-off error might propagate through the matrix multiplication, need some sort of error control
    for iStep=2:nSteps
      rho_t(:,:,:,iStep) = mmult(Q(:,:,:,iStep-1),...
                                   mmult(rho_t(:,:,:,iStep-1),...
                                         conj(permute(Q(:,:,:,iStep-1),[2,1,3,4])),'complex'),...
                                 'complex');                  
    end
    
    rho_t = rho_t(4:6,1:3,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'Oganesyan'  % see Ref [3]
    
    persistent cacheTensors
    
    rho_t = zeros(6,6,nTraj,nSteps);
    rho_t(1:3,4:6,:,1) = repmat(eye(3),1,1,nTraj,1);
    rho_t(4:6,1:3,:,1) = repmat(eye(3),1,1,nTraj,1);
%     rho_t = reshape(rho_t,[36,nTraj,nSteps]);

    if isempty(cacheTensors)
      % ISTOs in the lab frame and ISTs in the principal frame are
      % time-independent, so we only need to calculate them once
      
      SpinOps = cell(numel(Sys.Spins),3);
      for iSpin = 1:numel(Sys.Spins)
        SpinOps{iSpin,1} = sop(Sys.Spins,iSpin,1);
        SpinOps{iSpin,2} = sop(Sys.Spins,iSpin,2);
        SpinOps{iSpin,3} = sop(Sys.Spins,iSpin,3);
      end

      if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end  % TODO include frames in cardamom

      [T,F,~,~,~] = magint(Sys,SpinOps,CenterField,0,0);

      % sum over interactions
      F0 = F.F0(2)*2*pi;  % Hz -> rad s^-1, only keep isotropic HF interaction
      F2 = F.F2*2*pi;  % Hz -> rad s^-1

      T0 = T.T0{2};  % only keep isotropic HF interaction
      T2 = T.T2;
      
      cacheTensors.G0 = conj(F0)*T0;
      
      cacheTensors.G2 = zeros(6,6,5,5);
      
      for mp = 1:5
        for m = 1:5
          cacheTensors.G2(:,:,mp,m) = conj(F2(1,mp)).*T2{1,m} ...
                                      + conj(F2(2,mp)).*T2{2,m};
        end
      end
      
    end
    
    % calculate Wigner D-matrices
    [~,D2] = wigD(qTraj);
    
% -------------------------------------------------------------------------
    
%     H = cacheTensors.F0*cacheTensors.T0 ...
%         + squeeze(sum(cacheTensors.T2.*permute(conj(F2Traj),[4,5,1,2,3]),3));
    
%     U = zeros(6,6,nTraj,nSteps);
    
    for iStep=1:nSteps-1
      for iTraj=1:nTraj
%         U(:,:,iTraj,iStep) = expm(-1i*dt*H(:,:,iTraj,iStep));
        H = cacheTensors.G0;
        for mp = 1:5
          for m = 1:5
            H = H + D2(m,mp,iTraj,iStep).*cacheTensors.G2(:,:,mp,m);
          end
        end
%         U(:,:,iTraj,iStep) = expeig(1i*dt*H);
        U = expeig(1i*dt*H,1);
%         U = expm(1i*dt*H);
        rho_t(:,:,iTraj,iStep+1) = U*rho_t(:,:,iTraj,iStep)*U';
      end
    end
    
    % Only keep the m_S=-1/2 subspace part that contributes to 
    %   tr(S_{+}\rho(t))
    rho_t = rho_t(4:6,1:3,:,:);
    
%     Udag = conj(permute(U,[2,1,3,4]));
%         
%     for iStep=2:nSteps
%       rho_t(:,:,:,iStep) = mmult(U(:,:,:,iStep-1),...
%                                    mmult(rho_t(:,:,:,iStep-1),...
%                                          Udag(:,:,:,iStep-1),'complex'),...
%                                  'complex');                  
%     end
    
    % Superoperator code
%     % promote to Liouville space operator
%     L = tosuper(H,'c');
%     
% %     U = zeros(36,36,nTraj,nSteps);
%     
%     for iStep=2:nSteps
%       for iTraj=1:iTraj
% %         U = expm(-1i*dt*L(:,:,iTraj,iStep));
%         rho_t(:,iTraj,iStep) = expm(-1i*dt*L(:,:,iTraj,iStep))...
%                               *rho_t(:,iTraj,iStep-1);
%       end
%     end
%
%     rho_t = reshape(rho_t,[6,6,nTraj,nSteps]);

  otherwise
    error('Propagation method not recognized.')
end

end

% Helper functions
% -------------------------------------------------------------------------

function expmA = expeig(A,balancing)

if balancing==1
  [V,D] = eig(A);
elseif balancing==0  
  [V,D] = eig(A,'nobalance');
  V = V./sqrt(sum(V.*V,1));
end

expmA = V*diag(exp(diag(D)))*V';

end

function [D1,D2] = wigD(q)
% calculate Wigner D-matrices of specified rank from quaternions for 
% rotation of ISTOs

% nTraj = size(q,2);

A = q(1,:,:) - 1i*q(4,:,:);
B = -q(3,:,:) - 1i*q(2,:,:);
Ast = q(1,:,:) + 1i*q(4,:,:);
Bst = -q(3,:,:) + 1i*q(2,:,:);

Z = A.*Ast - B.*Bst;

% rank 1

% D1 = [            A.^2,       sqrt(2)*A.*B,           B.^2;
%        -sqrt(2)*A.*Bst,                  Z, sqrt(2)*Ast.*B;
%                 Bst.^2, -sqrt(2).*Ast.*Bst,         Ast.^2 ];      
D1(1,1,:,:) = A.^2;
D1(1,2,:,:) = sqrt(2)*A.*B;
D1(1,3,:,:) = B.^2;
D1(2,1,:,:) = -sqrt(2)*A.*Bst;
D1(2,2,:,:) = Z;
D1(2,3,:,:) = sqrt(2)*Ast.*B;
D1(3,1,:,:) = Bst.^2;
D1(3,2,:,:) = -sqrt(2).*Ast.*Bst;
D1(3,3,:,:) = Ast.^2;

% rank 2

% D2 = [                 A.^4,          2*A.^3.*B,     sqrt(6)*A.^2.*B.^2,         2*A.*B.^3,                 B.^4;
%                -2*A.^3.*Bst,      A.^2.*(2*Z-1),        sqrt(6)*A.*B.*Z,     B.^2.*(2*Z+1),          2*Ast.*B.^3;
%        sqrt(6)*A.^2.*Bst.^2, -sqrt(6)*A.*Bst.*Z,         1/2*(3*Z.^2-1), sqrt(6)*Ast.*B.*Z, sqrt(6)*Ast.^2.*B.^2;
%                -2*A.*Bst.^3,    Bst.^2.*(2*Z+1),   -sqrt(6)*Ast.*Bst.*Z,   Ast.^2.*(2*Z-1),          2*Ast.^3.*B;
%                      Bst.^4,     -2*Ast.*Bst.^3, sqrt(6)*Ast.^2.*Bst.^2,    -2*Ast.^3.*Bst,               Ast.^4 ];
D2(1,1,:,:) = A.^4;
D2(1,2,:,:) = 2*A.^3.*B;
D2(1,3,:,:) = sqrt(6)*A.^2.*B.^2;
D2(1,4,:,:) = 2*A.*B.^3;
D2(1,5,:,:) = B.^4;
D2(2,1,:,:) = -2*A.^3.*Bst;
D2(2,2,:,:) = A.^2.*(2*Z-1);
D2(2,3,:,:) = sqrt(6)*A.*B.*Z;
D2(2,4,:,:) = B.^2.*(2*Z+1);
D2(2,5,:,:) = 2*Ast.*B.^3;
D2(3,1,:,:) = sqrt(6)*A.^2.*Bst.^2;
D2(3,2,:,:) = -sqrt(6)*A.*Bst.*Z;
D2(3,3,:,:) = 1/2*(3*Z.^2-1);
D2(3,4,:,:) = sqrt(6)*Ast.*B.*Z;
D2(3,5,:,:) = sqrt(6)*Ast.^2.*B.^2;
D2(4,1,:,:) = -2*A.*Bst.^3;
D2(4,2,:,:) = Bst.^2.*(2*Z+1);
D2(4,3,:,:) = -sqrt(6)*Ast.*Bst.*Z;
D2(4,4,:,:) = Ast.^2.*(2*Z-1);
D2(4,5,:,:) = 2*Ast.^3.*B;
D2(5,1,:,:) = Bst.^4;
D2(5,2,:,:) = -2*Ast.*Bst.^3;
D2(5,3,:,:) = sqrt(6)*Ast.^2.*Bst.^2;
D2(5,4,:,:) = -2*Ast.^3.*Bst;
D2(5,5,:,:) = Ast.^4;


end