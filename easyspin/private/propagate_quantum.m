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
%
%   Opt: optional settings
%     Method         string
%                    'Sezer': propagate using the m_s=-1/2 subspace
%                    'DeSensi': propagate using an eigenvalue method
%                    'Oganesyan': propagate using correlation functions
%
%   MD:
%     RTraj          numeric, size = (3,3,nTraj,nSteps)
%                    externally provided rotation matrices
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

function rho_t = propagate_quantum(Sys, Par, Opt, MD, omega, CenterField)

% Preprocessing
% -------------------------------------------------------------------------

persistent cacheTensors

Method = Opt.Method;
Model = Par.Model;

if isfield(Par,'RTraj')
  RTraj = Par.RTraj;
  RTrajInv = permute(RTraj,[2,1,3,4]);
elseif isfield(Par,'qTraj')
  qTraj = Par.qTraj;
elseif isfield(MD,'RTraj')
  RTraj = MD.RTraj;
  RTrajInv = permute(RTraj,[2,1,3,4]);
else
  error('Par.RTraj, Par.qTraj, or MD.RTraj must be provided.')
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

    g_tr = sum(g);
    
    % Process MD simulation data
    % ---------------------------------------------------------------------
    
    if strcmp(Model,'Molecular Dynamics')
      % time step of MD simulation, MD.dt, is usually much smaller than
      % that of the propagation, Par.dt, so we need to average over windows
      % of size Par.dt/MD.dt
      
      if MD.nTraj > 1, error('Using the Sezer Method with multiple trajectories is not supported.'); end
      
      % size of averaging window
      nWindow = ceil(Par.dt/MD.dt);
      
      % size of MD trajectory after averaging
      M = floor(MD.nSteps/nWindow);

      GptensorAvg = zeros(3,3,M);
      AtensorAvg = zeros(3,3,M);
      
      % Perform rotations on g- and A-tensors
      Gptensor = matmult(RTraj, matmult(repmat(diag(g),1,1,MD.nTraj,MD.nSteps), ...
                                        RTrajInv))/gfree - g_tr/3/gfree;

      Atensor = matmult(RTraj, matmult(repmat(diag(A),1,1,MD.nTraj,MD.nSteps), ...
                                       RTrajInv))*1e6*2*pi;  % MHz (s^-1) -> Hz (rad s^-1)
    
      % average the interaction tensors over time windows
      for k = 1:M
        GptensorAvg(:,:,k) = mean(Gptensor(:,:,:,1+(k-1)*nWindow:k*nWindow),4);
        AtensorAvg(:,:,k) = mean(Atensor(:,:,:,1+(k-1)*nWindow:k*nWindow),4);
      end
      
      % process single long trajectory into multiple short trajectories
      lag = 1;
      if Par.nSteps<M
        nSteps = Par.nSteps;
        nTraj = floor((M-nSteps)/lag) + 1;
      else
        nSteps = M;
        nTraj = 1;
      end
      
      Gptensor = zeros(3,3,nTraj,nSteps);
      Atensor = zeros(3,3,nTraj,nSteps);
      
      for k = 1:nTraj
        idx = (1:nSteps) + (k-1)*lag;
        Gptensor(:,:,k,:) = GptensorAvg(:,:,idx);
        Atensor(:,:,k,:) = AtensorAvg(:,:,idx);
      end
      
      % frame of MD coordinate system is lab frame for solution simulations
      
    else

      % Perform rotations on g- and A-tensors
      Gptensor = matmult(RTraj, matmult(repmat(diag(g),1,1,nTraj,nSteps), ...
                                            RTrajInv))/gfree - g_tr/3/gfree;

      Atensor = matmult(RTraj, matmult(repmat(diag(A),1,1,nTraj,nSteps), ...
                                           RTrajInv))*1e6*2*pi;  % MHz (s^-1) -> Hz (rad s^-1)
      
    end
    
    rho_t = zeros(3,3,nTraj,nSteps);
    rho_t(:,:,:,1) = repmat(eye(3),[1,1,nTraj]);
    
    % Prepare propagators
    % ---------------------------------------------------------------------
    
    Gp_zz = Gptensor(3,3,:,:);

    % norm of expression in Eq. 24 in [1]
    a = sqrt(Atensor(1,3,:,:).*Atensor(1,3,:,:) ...
           + Atensor(2,3,:,:).*Atensor(2,3,:,:) ...
           + Atensor(3,3,:,:).*Atensor(3,3,:,:));

    % rotation angle and unit vector parallel to axis of rotation
    % refer to paragraph below Eq. 37 in [1]
%     theta = Gamma*dt*0.5*squeeze(a);
    theta = dt*0.5*squeeze(a);
    nx = squeeze(Atensor(1,3,:,:)./a);
    ny = squeeze(Atensor(2,3,:,:)./a);
    nz = squeeze(Atensor(3,3,:,:)./a);

    % Eqs. A1-A2 in [1]
    ct = cos(theta) - 1;
    st = -sin(theta);

    % matrix exponential of hyperfine part
    % Eqs. A1-A2 are used to construct Eq. 37 in [1]
    expadotI = zeros(3,3,size(Gp_zz,3),size(Gp_zz,4));
    expadotI(1,1,:,:) = 1 + ct.*(nz.*nz + 0.5*(nx.*nx + ny.*ny)) ...
                         + 1i*st.*nz;
    expadotI(1,2,:,:) = sqrt(0.5)*(st.*ny + ct.*nz.*nx) ...
                         + 1i*sqrt(0.5)*(st.*nx - ct.*nz.*ny);
    expadotI(1,3,:,:) = 0.5*ct.*(nx.*nx - ny.*ny) ... 
                         - 1i*ct.*nx.*ny;
    expadotI(2,1,:,:) = sqrt(0.5)*(-st.*ny + ct.*nz.*nx) ...
                         + 1i*sqrt(0.5)*(st.*nx + ct.*nz.*ny);
    expadotI(2,2,:,:) = 1 + ct.*(nx.*nx + ny.*ny);
    expadotI(2,3,:,:) = sqrt(0.5)*(st.*ny - ct.*nz.*nx) ...
                         + 1i*sqrt(0.5)*(st.*nx + ct.*nz.*ny);
    expadotI(3,1,:,:) = 0.5*ct.*(nx.*nx - ny.*ny) ...
                         + 1i*ct.*nx.*ny;
    expadotI(3,2,:,:) = sqrt(0.5)*(-st.*ny - ct.*nz.*nx) ...
                         + 1i*sqrt(0.5)*(st.*nx - ct.*nz.*ny);
    expadotI(3,3,:,:) = 1 + ct.*(nz.*nz + 0.5*(nx.*nx + ny.*ny))  ...
                         - 1i*st.*nz;


    % Calculate propagator, Eq. 35 in [1]
    U = bsxfun(@times, exp(-1i*dt*0.5*omega*Gp_zz), expadotI);

    
    % Propagate density matrix
    % ---------------------------------------------------------------------
    
%     if strcmp(Model,'Molecular Dynamics')
%       % only one trajectory, so the "*" operator can be used
%       for iStep=2:M
%         rho_t(:,:,:,iStep) = U(:,:,:,iStep-1)*rho_t(:,:,:,iStep-1)*U(:,:,:,iStep-1);
%       end
%     else
      % there are multiple trajectories, so we need "mmult"
    if nTraj>1
      for iStep=2:nSteps
        rho_t(:,:,:,iStep) = mmult(U(:,:,:,iStep-1),...
                                   mmult(rho_t(:,:,:,iStep-1), U(:,:,:,iStep-1),'complex'),...
                                   'complex');

        %  Trace of a product of matrices is the sum of entry-wise products
        %      rho_t(:,:,:,iStep) = U2(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
      end
    else
      for iStep=2:nSteps
        rho_t(:,:,1,iStep) = U(:,:,1,iStep-1)*rho_t(:,:,1,iStep-1)*U(:,:,1,iStep-1);

        %  Trace of a product of matrices is the sum of entry-wise products
        %      rho_t(:,:,:,iStep) = U2(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
      end
    end
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'DeSensi'  % see Ref [2]

    rho_t = zeros(6,6,nTraj,nSteps);
    rho_t(1:3,4:6,:,1) = repmat(eye(3),1,1,nTraj,1);
    rho_t(4:6,1:3,:,1) = repmat(eye(3),1,1,nTraj,1);

    g_tr = sum(g);
    
    B0 = CenterField/1e3;  % mT -> T

    % Perform rotations on g- and A-tensors
    gtensor = matmult(RTraj, matmult(repmat(diag(g),1,1,nTraj,nSteps), ...
                                      RTrajInv));

    Atensor = matmult(RTraj, matmult(repmat(diag(A),1,1,nTraj,nSteps), ...
                                      RTrajInv))*1e6*2*pi;  % MHz (s^-1) -> rad s^-1
                                    
    % Prepare propagators
    % ---------------------------------------------------------------------

    g_zz = gtensor(3,3,:,:);

    A_xz = Atensor(1,3,:,:);
    A_yz = Atensor(2,3,:,:);
    A_zz = Atensor(3,3,:,:);

    % notation adapted from Eq. 24-27 in [2]
    omegaeff = -bmagn/(planck/2/pi)*B0*(g_zz - g_tr/3); 
    
    ma = 0.5;
    mb = -0.5;
    
    mp = 1;
    m0 = 0;
    mm = -1;

    ca = A_zz;% - 4*ma*omegaN;
    cb = A_zz;% - 4*mb*omegaN;

    ella = sqrt( (A_xz).^2 + (A_yz).^2 + ca.^2 );
    ellb = sqrt( (A_xz).^2 + (A_yz).^2 + cb.^2 );

    % Matrix of eigenvalues, adapted from Eqs. 24-27 in Ref. [2]
    expLambda = zeros(6,6,nTraj,nSteps);
    expLambda(1,1,:,:) = exp(-1i*dt*ma*(omegaeff + mp*ella));
    expLambda(2,2,:,:) = exp(-1i*dt*ma*(omegaeff + m0*ella));
    expLambda(3,3,:,:) = exp(-1i*dt*ma*(omegaeff + mm*ella));
    expLambda(4,4,:,:) = exp(-1i*dt*mb*(omegaeff + mp*ellb));
    expLambda(5,5,:,:) = exp(-1i*dt*mb*(omegaeff + m0*ellb));
    expLambda(6,6,:,:) = exp(-1i*dt*mb*(omegaeff + mm*ellb));

    % Matrix of eigenvectors, adapted from Eq. 28 in Ref. [2]
    % NOTE: b (and bst) are not defined in the text, and the expressions in
    % Eq. 28 are not correct; as was below, the following factors need to
    % be replaced as follows:
    %         1/4 --> 1
    %   1/sqrt(2) --> sqrt(2)
    
    b   = (A_xz + 1i*A_yz);
    bst = (A_xz - 1i*A_yz);
    
    V  = [ (ca+ella).^2./b.^2,                     -bst./b,      (ca-ella).^2./b.^2, zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps);
           sqrt(2)*(ca+ella)./b,             sqrt(2)*ca./b,    sqrt(2)*(ca-ella)./b, zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps);
           ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps);
          zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps),      (cb+ellb).^2./b.^2,                 -bst./b,      (cb-ellb).^2./b.^2;
          zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps),    sqrt(2)*(cb+ellb)./b,           sqrt(2)*cb./b,    sqrt(2)*(cb-ellb)./b;
          zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps), zeros(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),  ones(1,1,nTraj,nSteps),   ones(1,1,nTraj,nSteps) ];

    % Normalize eigenvectors
    V = bsxfun(@rdivide,V,sqrt(sum(V.*conj(V),1)));
    
    % Calculate propagator
    U = mmult(V, mmult(expLambda, conj(permute(V,[2,1,3,4])), 'complex'), 'complex');
    
    Udag = conj(permute(U,[2,1,3,4]));
    
    % Propagate density matrix
    % ---------------------------------------------------------------------
    
% FIXME round-off error might propagate through the matrix multiplication, need some sort of error control
    for iStep=2:nSteps
      rho_t(:,:,:,iStep) = mmult(U(:,:,:,iStep-1),...
                                   mmult(rho_t(:,:,:,iStep-1),...
                                         Udag(:,:,:,iStep-1),'complex'),...
                                 'complex');                  
    end
    
    rho_t = rho_t(4:6,1:3,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'Oganesyan'  % see Ref [3]
    
    rho_t = zeros(6,6,nTraj,nSteps);
    rho_t(1:3,4:6,:,1) = repmat(eye(3),1,1,nTraj,1);
    rho_t(4:6,1:3,:,1) = repmat(eye(3),1,1,nTraj,1);
%     rho_t = zeros(6,6,1,nSteps);
%     rho_t(1:3,4:6,1,1) = repmat(eye(3),1,1,1,1);
%     rho_t(4:6,1:3,1,1) = repmat(eye(3),1,1,1,1);

    % Calculate and store rotational basis operators
    % ---------------------------------------------------------------------
    if isempty(cacheTensors)
      % ISTOs in the lab frame and IST components in the principal frame 
      % are time-independent, so we only need to calculate them once
      
      SpinOps = cell(numel(Sys.Spins),3);
%       for iSpin = 1:numel(Sys.Spins)
%         SpinOps{iSpin,1} = sop(Sys.Spins,iSpin,1);
%         SpinOps{iSpin,2} = sop(Sys.Spins,iSpin,2);
%         SpinOps{iSpin,3} = sop(Sys.Spins,iSpin,3);
%       end
      SpinOps{1,1} = zeros(6);
      SpinOps{1,2} = zeros(6);
      SpinOps{1,3} = sop(Sys.Spins,1,3);
      
      SpinOps{2,1} = sop(Sys.Spins,2,1);
      SpinOps{2,2} = sop(Sys.Spins,2,2);
      SpinOps{2,3} = sop(Sys.Spins,2,3);
      
      if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end  % TODO include frames in cardamom
      
      [T,F,~,~,~] = magint(Sys,SpinOps,CenterField,0,0);

      F0 = F.F0(2)*2*pi;  % Hz -> rad s^-1, only keep isotropic HF interaction
      F2 = F.F2*2*pi;  % Hz -> rad s^-1

      T0 = T.T0{2};  % only keep isotropic HF interaction
      T2 = T.T2;
      
      % zeroth rank
      cacheTensors.Q0 = conj(F0)*T0;
      
      cacheTensors.Q2 = cell(5,5);
      
      % create the 25 second-rank RBOs
      for mp = 1:5
        for m = 1:5
          cacheTensors.Q2{mp,m} = conj(F2(1,mp))*T2{1,m} ...
                                      + conj(F2(2,mp))*T2{2,m};
        end
      end
      
    end
    
    % Prepare propagators
    % ---------------------------------------------------------------------
    
    % calculate Wigner D-matrices from the quaternion trajectories
    [~, D2] = wigD(qTraj);
    

    H = zeros(6,6,nTraj,nSteps);
    
    % zeroth rank term (HF only)
%     H(:,:,:,:) = repmat(cacheTensors.Q0,[1,1,nTraj,nSteps]);
    
    % calculate Hamiltonians
    for iStep=1:nSteps
      for iTraj=1:nTraj
        % zeroth rank term (HF only)
        H(:,:,iTraj,iStep) = cacheTensors.Q0;
        
        % rotate second rank terms and add to Hamiltonian
        for mp = 1:5
          for m = 1:5
            H(:,:,iTraj,iStep) = H(:,:,iTraj,iStep) + D2(m,mp,iTraj,iStep)*cacheTensors.Q2{mp,m};
          end
        end
      end
    end
%     Hmean = mean(H,4);
%     deltaH = H - Hmean;

    U = zeros(6,6,nTraj,nSteps);

    % calculate propagators
    for iStep=1:nSteps
      for iTraj=1:nTraj
%         U(:,:,iTraj,iStep) = expm(-1i*dt*H);
        U(:,:,iTraj,iStep) = expeig(1i*dt*H(:,:,iTraj,iStep));
      end
    end
    
%     U = U - mean(mean(U,3),4);
    
    Udag = conj(permute(U,[2,1,3,4]));
        
    % Propagate density matrix
    % ---------------------------------------------------------------------
    for iStep=2:nSteps
      rho_t(:,:,:,iStep) = mmult(U(:,:,:,iStep-1),...
                                 mmult(rho_t(:,:,:,iStep-1),...
                                       Udag(:,:,:,iStep-1),'complex'),...
                                 'complex');                  
    end
%     for iStep=2:nSteps
%       rho_t(:,:,:,iStep) = U(:,:,1,iStep-1)*rho_t(:,:,1,iStep-1)...
%                                            *U(:,:,1,iStep-1)';                  
%     end
    
    % Only keep the m_S=-1/2 subspace part that contributes to 
    %   tr(S_{+}\rho(t))
    rho_t = rho_t(4:6,1:3,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise
    error('Propagation method not recognized.')
end


end

% Helper functions
% -------------------------------------------------------------------------

function C = expeig(A)

[V,D] = eig(A);

C = V*diag(exp(diag(D)))*V';

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
D1 = zeros(3,3,size(q,2),size(q,3));
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
D2 = zeros(5,5,size(q,2),size(q,3));
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