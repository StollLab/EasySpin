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
% [2] Oganesyan, Phys. Chem. Chem. Phys. 13, 4724 (2011)
%      http://dx.doi.org/10.1039/c0cp01068e

function rho_t = propagate_quantum(Sys, Par, Opt, MD, omega, CenterField)

% Preprocessing
% -------------------------------------------------------------------------

persistent cacheTensors
persistent qTraj
persistent fullSteps
persistent K2

if ~isfield(Opt, 'Liouville')
  Opt.Liouville = 0;
end

Liouville = Opt.Liouville;
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
if ~isfield(Par,'truncated'), Par.truncated = 1; end

truncated = Par.truncated;

dt = Par.dt;
nTraj = Par.nTraj;
if strcmp(Model, 'Molecular Dynamics')
  nSteps = Par.nSteps;
else
  if strcmp(Method, 'Nitroxide')
    nSteps = size(RTraj, 4);
  elseif strcmp(Method, 'ISTOs')
    nSteps = size(qTraj, 3);
  end
end
t = linspace(0, dt*nSteps, nSteps);

if ~isequal(size(g),[1,3]) || ~isequal(size(A),[1,3])
  error('g and A tensor values must be 3-vectors.')
end

% Set up for MD time averaging

if strcmp(Model,'Molecular Dynamics')
  % time step of MD simulation, MD.dt, is usually much smaller than
  % that of the propagation, Par.dt, so determine the size of the averaging
  % window here (Par.dt/MD.dt)
  
  
  if MD.nTraj > 1, error('Using multiple MD trajectories is not supported.'); end
  
  if any(strcmp(Method,{'Nitroxide','ISTOs'}))
    % size of averaging window
    nWindow = ceil(Par.dt/MD.dt);

    % size of MD trajectory after averaging
    M = floor(MD.nSteps/nWindow);

    % process single long trajectory into multiple short trajectories
%     lag = ceil(Par.dt/2e-9);  % use 2 ns lag between windows
      lag = 2;
%     lag = ceil(Par.dt/1e-9);
    if Par.nSteps<M
      nSteps = Par.nSteps;
      nTraj = floor((M-nSteps)/lag) + 1;
    else
      
      nSteps = M;
      nTraj = 1;
    end

  elseif strcmp(Method,'Resampling')
    % Steinhoff method does not require tensor averaging, so set nTraj to
    % user input
    nTraj = Par.nTraj;
  end
  
  if isfield(MD,'GlobalDiff')
    % generate global isotropic diffusion
    
    if isfield(Sys,'Coefs'), Sys = rmfield(Sys,'Coefs'); end
    if isfield(Sys,'PseudoPotFun'), Sys = rmfield(Sys, 'PseudoPotFun'); end
    
    Sys.Diff = MD.GlobalDiff;
    Par.nTraj = nTraj;
    Par.nSteps = nSteps;
    [t, RTrajGlobal, qTrajGlobal] = stochtraj(Sys,Par);
    RTrajGlobalInv = permute(RTrajGlobal,[2,1,3,4]);
  end
  
end

% Gamma = gfree*bmagn/(planck/2/pi);  % rad s^-1 T^-1
% omegaN = 19.331e6*B;  % gyromagnetic ratio for 14N: 
                      % 19.331x10^6 rad s^-1 T^-1

% Simulation
% -------------------------------------------------------------------------

switch Method
  case 'Nitroxide'  % see Ref [1]

    g_tr = sum(g);
    
    % Process MD simulation data
    % ---------------------------------------------------------------------
    
    % Perform rotations on g- and A-tensors
    gTensor = tensortraj(g,RTraj,RTrajInv);

    ATensor = tensortraj(A,RTraj,RTrajInv)*1e6*2*pi; % MHz (s^-1) -> Hz (rad s^-1)
    
    if strcmp(Model,'Molecular Dynamics')

      gTensorAvg = zeros(3,3,M);
      ATensorAvg = zeros(3,3,M);
    
      % average the interaction tensors over time windows
      idx = 1:nWindow;
      for k = 1:M
        gTensorAvg(:,:,k) = mean(gTensor(:,:,:,idx),4);
        ATensorAvg(:,:,k) = mean(ATensor(:,:,:,idx),4);
        idx = idx + nWindow;
      end
      
      gTensor = zeros(3,3,nTraj,nSteps);
      ATensor = zeros(3,3,nTraj,nSteps);
      
      for k = 1:nTraj
        idx = (1:nSteps) + (k-1)*lag;
        gTensor(:,:,k,:) = gTensorAvg(:,:,idx);
        ATensor(:,:,k,:) = ATensorAvg(:,:,idx);
      end
      
      if isfield(MD,'GlobalDiff')
        gTensor = matmult(RTrajGlobal, matmult(gTensor, RTrajGlobalInv));
        ATensor = matmult(RTrajGlobal, matmult(ATensor, RTrajGlobalInv));
      end
        
    end
    
    GpTensor = gTensor/gfree - g_tr/3/gfree;
        
    rho_t = zeros(3,3,nTraj,nSteps);
    rho_t(:,:,:,1) = repmat(eye(3),[1,1,nTraj]);
    
    % Prepare propagators
    % ---------------------------------------------------------------------
    
    Gp_zz = GpTensor(3,3,:,:);

    % norm of expression in Eq. 24 in [1]
    a = sqrt(ATensor(1,3,:,:).*ATensor(1,3,:,:) ...
           + ATensor(2,3,:,:).*ATensor(2,3,:,:) ...
           + ATensor(3,3,:,:).*ATensor(3,3,:,:));

    % rotation angle and unit vector parallel to axis of rotation
    % refer to paragraph below Eq. 37 in [1]
%     theta = Gamma*dt*0.5*squeeze(a);
    theta = dt*0.5*squeeze(a);
    nx = squeeze(ATensor(1,3,:,:)./a);
    ny = squeeze(ATensor(2,3,:,:)./a);
    nz = squeeze(ATensor(3,3,:,:)./a);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'ISTOs'  % see Ref [2]

    % Calculate and store rotational basis operators
    % ---------------------------------------------------------------------
    if isempty(cacheTensors)
      % ISTOs in the lab frame and IST components in the principal frame 
      % are time-independent, so we only need to calculate them once
      
      
      % set up spin operators in |S,m_S>|I,m_I> product space

      % electron spin operators
%       SpinOps{1,1} = sop(Sys.Spins,1,1);
%       SpinOps{1,2} = sop(Sys.Spins,1,2);
      SpinOps{1,1} = zeros(6);  % S_x and S_y are zero operators in HF limit
      SpinOps{1,2} = zeros(6);
      SpinOps{1,3} = sop(Sys.Spins,1,3);
      
      % nuclear spin operators
      SpinOps{2,1} = sop(Sys.Spins,2,1);
      SpinOps{2,2} = sop(Sys.Spins,2,2);
      SpinOps{2,3} = sop(Sys.Spins,2,3);
      
      if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end  % TODO include frames in cardamom
      
      [T,F,~,~,~] = magint(Sys,SpinOps,CenterField,0,0);

%       F0 = F.F0*2*pi;
      F0 = F.F0(2)*2*pi;  % Hz -> rad s^-1, only keep isotropic HF interaction
      F2 = F.F2*2*pi;  % Hz -> rad s^-1

%       T0 = T.T0;
      T0 = T.T0{2};  % only keep isotropic HF interaction
      T2 = T.T2;
      
      % zeroth rank
%       cacheTensors.Q0 = conj(F0(1))*T0{1} + conj(F0(2))*T0{2};

%       if Liouville
%         cacheTensors.Q0 = tosuper(conj(F0)*T0,'c');
%       else
        cacheTensors.Q0 = conj(F0)*T0;
%       end
      
      cacheTensors.Q2 = cell(5,5);
      
      % create the 25 second-rank RBOs
      for mp = 1:5
        for m = 1:5
%           if Liouville
%             cacheTensors.Q2{mp,m} = tosuper(conj(F2(1,mp))*T2{1,m} ...
%                                 + conj(F2(2,mp))*T2{2,m},'c');
%           else
            cacheTensors.Q2{mp,m} = conj(F2(1,mp))*T2{1,m} ...
                                        + conj(F2(2,mp))*T2{2,m};
%           end
        end
      end
      
    end
    
    % Prepare explicit propagators
    % ---------------------------------------------------------------------
    
    % MD trajectories yield rotation matrices, so convert to quaternions
    if strcmp(Model,'Molecular Dynamics')&&isempty(qTraj)
      qTraj = rotmat2quat(RTraj);
    end
    
    % calculate Wigner D-matrices from the quaternion trajectories
    D2 = wigD(qTraj);
    
    if truncated
      if isempty(fullSteps)&&isempty(K2)
        % set time step for integration of correlation functions
        if strcmp(Model,'Molecular Dynamics')
          % for MD trajectories, use MD timestep
          dtau = MD.dt;
        else
          % for stochastic dynamics, use propagation time step
          dtau = Par.dt;
        end

        % use size(D2,3) since the full MD trajectories need to be
        % integrated, not the windowed trajectories
        acorrD200 = zeros(size(D2,3), size(D2,4));
        for iTraj=1:size(D2,3)
          % normalized autocorrelation functions are needed here
          acorrD200(iTraj,:) = autocorrfft(squeeze(D2(3,3,iTraj,:)).', 0, 1, 1);
        end

        acorrD200 = mean(acorrD200, 1);
        time = linspace(0, size(D2,4)*dtau, size(D2,4));

        % calculate correlation time
        tauc = max(cumtrapz(time, acorrD200));  % FIXME find a robust way to calculate this

        % set the number of time steps to use full propagator before 
        % correlation functions relax
        fullSteps = ceil(truncated*tauc/dt);

        % calculate correlation functions for approximate propagator
        acorr = zeros(size(D2,3), size(D2,4));
        K2 = zeros(5,5,nTraj);
        for mp=1:5
          for m=1:5
            % non-normalized autocorrelation functions are needed here
            for iTraj=1:size(D2,3)
              acorr(iTraj,:) = autocorrfft(squeeze(D2(mp,m,iTraj,:)).',0,0,1);
            end
            K2(mp,m) = max(max(cumtrapz(time, mean(acorr,1))));  % TODO implement cross-correlation functions
          end
        end

        % if correlation time is longer than user input, just use full 
        % propagation scheme the entire time
        if fullSteps>nSteps
          fullSteps = nSteps; 
          truncated = 0;
        end
      end
      
    else
      % no approximation, use full propagator for all time steps
      fullSteps = nSteps;
      
    end
    
    if strcmp(Model,'Molecular Dynamics')

      D2Avg = zeros(5,5,M);
    
      % average the interaction tensors over time windows
      idx = 1:nWindow;
      for k = 1:M
        D2Avg(:,:,k) = mean(D2(:,:,:,idx),4);
        idx = idx + nWindow;
      end
      
      D2 = zeros(5,5,nTraj,nSteps);
      
      for k = 1:nTraj
        idx = (1:nSteps) + (k-1)*lag;
        D2(:,:,k,:) = D2Avg(:,:,idx);
      end
      
      if isfield(MD,'GlobalDiff')
        D2Global = wigD(qTrajGlobal);
        
        % combine Wigner D-matrices of global and averaged local motion
        D2 = mmult(D2Global, D2, 'complex');
      end
    end
    
    rho_t = zeros(6,6,nTraj,nSteps);
    rho_t(1:3,4:6,:,1) = repmat(eye(3),1,1,nTraj,1);
    rho_t(4:6,1:3,:,1) = repmat(eye(3),1,1,nTraj,1);
    
    % calculate Hamiltonians and  propagators
%     if Liouville
%       rho_t = reshape(rho_t,[36,nTraj,nSteps]);
%       U = zeros(36,36,nTraj,nSteps);
%     else
      U = zeros(6,6,nTraj,nSteps);
%     end

    if Liouville
      rho_t = reshape(rho_t,[36,nTraj,nSteps]);
    end
    
%     if truncated
%       AutoCorrFFT = autocorrfft(squeeze(D2(3,3,:,:)), 0, 1, 1);
% 
%       N = round(nSteps/2);
%       M = round(N/2);
% 
%       AutoCorrFFT = mean(AutoCorrFFT, 1).';
% 
%       % the ACF will not always reach zero for a small number of trajectories,
%       % so subtract the offset, which does not change the correlation time
%       AutoCorrFFT = AutoCorrFFT - mean(AutoCorrFFT(M:3*M));
%       
% %       AutoCorrFFTavg = AutoCorrFFTavg + AutoCorrFFT;
% 
%       % calculate correlation time
%       tau = trapz(t(1:N), AutoCorrFFT(1:N));  % FIXME find a robust way to calculate this
%       if tau>tauc, tauc = tau; end
% 
%       fullSteps = 10*ceil(tauc/dt);  % build full propagator until 
%                                      % correlation functions have relaxed
%                                      
%       % if correlation time is longer than user input, just use full 
%       % propagation scheme the entire time
%       if fullSteps>nSteps
%         fullSteps = nSteps; 
%         truncated = 0;
%       end
%       
%     else
%       fullSteps = nSteps;  % build full propagator for all time steps
%       
%     end
    
    for iStep=1:fullSteps
      for iTraj=1:nTraj
        % zeroth rank term (HF only)
%         H(:,:,iTraj,iStep) = cacheTensors.Q0;
        H = cacheTensors.Q0;
%         Hfull(:,:,iTraj,iStep) = H;
        
        % rotate second rank terms and add to Hamiltonian
        for mp = 1:5
          for m = 1:5
%             Hfull(:,:,iTraj,iStep) = Hfull(:,:,iTraj,iStep) + D2(m,mp,iTraj,iStep)*cacheTensors.Q2{mp,m};
            H = H + D2(m,mp,iTraj,iStep)*cacheTensors.Q2{mp,m};
          end
        end
        U(:,:,iTraj,iStep) = expeig(1i*dt*H);  % TODO speed this up!
%         U(:,:,iTraj,iStep) = expm_fast1(1i*dt*H);  % TODO speed this up!
      end
    end
    
    Udag = conj(permute(U,[2,1,3,4]));
    
    if Liouville
      U = tosuperLR(U, Udag);
    end
    
    % Propagate density matrix
    % ---------------------------------------------------------------------
    
    if truncated
      % Prepare equilibrium propagators
      % -------------------------------------------------------------------
      
%       if Liouville
%         Ueq = zeros(36,36,nTraj);
%       else
        Ueq = zeros(6,6,nTraj);
%         Ueqdag = zeros(6,6,nTraj);
%       end
      
      for iTraj=1:nTraj
        % zeroth rank term (HF only)
        Heq = cacheTensors.Q0;
        
        % rotate second rank terms and add to Hamiltonian
        for mp = 1:5
          for m = 1:5
            Heq = Heq + cacheTensors.Q2{mp,m}*mean(D2(m,mp,iTraj,:),4) ...
                      + 1i*cacheTensors.Q2{mp,m}*(cacheTensors.Q2{mp,m})'*K2(mp,m);
          end
        end
%         Ueq(:,:,iTraj) = expeig(1i*dt*Heq);
        Ueq(:,:,iTraj) = expm_fast1(1i*dt*Heq);
      end
      
      Ueqdag = conj(permute(Ueq,[2,1,3]));
      
      if Liouville
        Ueq = tosuperLR(Ueq,Ueqdag);
      end
      
      % full propagation until correlation functions have relaxed
      if Liouville
        for iStep=2:fullSteps
          for iTraj=1:nTraj
            rho_t(:,iTraj,iStep) = U(:,:,iTraj,iStep-1)*rho_t(:,iTraj,iStep-1);  
          end
        end

        % equilibrium propagation
        for iStep=fullSteps+1:nSteps
          for iTraj=1:nTraj
            rho_t(:,iTraj,iStep) = Ueq(:,:,iTraj)*rho_t(:,iTraj,iStep-1);  
          end                 
        end
        
      else
        for iStep=2:fullSteps
          rho_t(:,:,:,iStep) = mmult(U(:,:,:,iStep-1),...
                                     mmult(rho_t(:,:,:,iStep-1),...
                                           Udag(:,:,:,iStep-1),'complex'),...
                                     'complex');                  
        end

        % equilibrium propagation
        for iStep=fullSteps+1:nSteps
          rho_t(:,:,:,iStep) = mmult(Ueq,...
                                     mmult(rho_t(:,:,:,iStep-1),...
                                           Ueqdag,'complex'),...
                                     'complex');                  
        end
      end
      
    else
    % truncated trajectory behavior not being used, so use full propagation
    % scheme for all time steps
      if Liouville
        for iStep=2:nSteps
          for iTraj=1:nTraj
            rho_t(:,iTraj,iStep) = U(:,:,iTraj,iStep-1)*rho_t(:,iTraj,iStep-1);  
          end
        end
        
      else
        for iStep=2:nSteps
          rho_t(:,:,:,iStep) = mmult(U(:,:,:,iStep-1),...
                                     mmult(rho_t(:,:,:,iStep-1),...
                                           Udag(:,:,:,iStep-1),'complex'),...
                                     'complex');                  
        end
      end
    end

    if Liouville
      rho_t = reshape(rho_t,[6,6,nTraj,nSteps]);
    end
    
    % Only keep the m_S=-1/2 subspace part that contributes to 
    %   tr(S_{+}\rho(t))
    rho_t = rho_t(4:6,1:3,:,:);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'Resampling'
    
    g_tr = sum(g);
    
    % Process MD simulation data
    % ---------------------------------------------------------------------
    
    % Perform rotations on g- and A-tensors
    gTensor = tensortraj(g,RTraj,RTrajInv);

    ATensor = tensortraj(A,RTraj,RTrajInv)*1e6*2*pi; % MHz (s^-1) -> Hz (rad s^-1)
    
    if isfield(MD,'GlobalDiff')
      gTensor = mmult(RTrajGlobal, mmult(gTensor, RTrajGlobalInv, 'real'), 'real');
      ATensor = mmult(RTrajGlobal, mmult(ATensor, RTrajGlobalInv, 'real'), 'real');
    end
    
    GpTensor = gTensor/gfree - g_tr/3/gfree;
        
    rho_t = zeros(3,3,Par.nTraj,nSteps);
    rho_t(:,:,:,1) = repmat(eye(3),[1,1,Par.nTraj]);
    
    % Prepare propagators
    % ---------------------------------------------------------------------
    
    Gp_zz = GpTensor(3,3,:,:);

    % norm of expression in Eq. 24 in [1]
    a = sqrt(ATensor(1,3,:,:).*ATensor(1,3,:,:) ...
           + ATensor(2,3,:,:).*ATensor(2,3,:,:) ...
           + ATensor(3,3,:,:).*ATensor(3,3,:,:));

    % rotation angle and unit vector parallel to axis of rotation
    % refer to paragraph below Eq. 37 in [1]
%     theta = Gamma*dt*0.5*squeeze(a);
    theta = dt*0.5*squeeze(a);
    nx = squeeze(ATensor(1,3,:,:)./a);
    ny = squeeze(ATensor(2,3,:,:)./a);
    nz = squeeze(ATensor(3,3,:,:)./a);

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

function F = expm_fast1(A)
%EXPM  Matrix exponential.
%   EXPM(A) is the matrix exponential of A and is computed using
%   a scaling and squaring algorithm with a Pade approximation.
%
%   Although it is not computed this way, if A has a full set
%   of eigenvectors V with corresponding eigenvalues D then
%   [V,D] = EIG(A) and EXPM(A) = V*diag(exp(diag(D)))/V.
%
%   EXP(A) computes the exponential of A element-by-element.
%
%   See also LOGM, SQRTM, FUNM.

%   References:
%   N. J. Higham, The scaling and squaring method for the matrix
%      exponential revisited. SIAM J. Matrix Anal. Appl., 26(4), (2005),
%      pp. 1179-1193.
%   A. H. Al-Mohy and N. J. Higham, A new scaling and squaring algorithm
%      for the matrix exponential, SIAM J. Matrix Anal. Appl., 31(3),
%      (2009), pp. 970-989.
%
%   Nicholas J. Higham and Samuel D. Relton
%   Copyright 1984-2015 The MathWorks, Inc.

T = A;

if isdiag(T) % Check if T is diagonal.
  d = diag(T);
  F = diag(exp(full(d)));
  return;
end

% Compute exponential
% Get scaling and Pade parameters.
[s, m, T2, T4, T6] = expm_params(T);

% Rescale the powers of T appropriately.
if s ~= 0
  T = T/(2.^s);
  T2 = T2/2^(s*2);
  T4 = T4/2^(s*4);
  T6 = T6/2^(s*6);
end

% Evaluate the Pade approximant.
switch m
  case 3
    c = [120, 60, 12, 1];
  case 5
    c = [30240, 15120, 3360, 420, 30, 1];
  case 7
    c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
  case 9
    c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
      2162160, 110880, 3960, 90, 1];
  case 13
    c = [64764752532480000, 32382376266240000, 7771770303897600, ...
      1187353796428800,  129060195264000,   10559470521600, ...
      670442572800,      33522128640,       1323241920,...
      40840800,          960960,            16380,  182,  1];
end
I = eye(length(T));
switch m
  case {3, 5, 7, 9}
    Tpowers = {[],T2,[],T4,[],T6};
    strt = length(Tpowers) + 2;
    for k = strt:2:m-1
      Tpowers{k} = Tpowers{k-2}*Tpowers{2};
    end
    U = c(2)*I;
    V = c(1)*I;
    for j = m:-2:3
      U = U + c(j+1)*Tpowers{j-1};
      V = V + c(j)*Tpowers{j-1};
    end
    U = T*U;
  case 13
    U = T * (T6*(c(14)*T6 + c(12)*T4 + c(10)*T2) + c(8)*T6 + c(6)*T4 + c(4)*T2 + c(2)*I);
    V = T6*(c(13)*T6 + c(11)*T4 + c(9)*T2) + c(7)*T6 + c(5)*T4 + c(3)*T2 + c(1)*I;
end
F = (V-U)\(2*U) + I;  %F = (-U+V)\(U+V);

% Squaring phase.
for k = 1:s
  F = F*F;
end
end


function t = ell(T, coeff, m_val)
%ell Function needed to compute optimal parameters.
scaledT = coeff.^(1/(2*m_val+1)) .* abs(T);
alpha = normAm(scaledT,2*m_val+1)/norm(T,1);
t = max(ceil(log2(2*alpha/eps(class(alpha)))/(2*m_val)),0);
end

function [s, m, T2, T4, T6] = expm_params(T)
%expm_params Obtain scaling parameter and order of the Pade approximant.
% Coefficients of backwards error function.
coeff = [1/100800, 1/10059033600, 1/4487938430976000,...
  1/5914384781877411840000, 1/113250775606021113483283660800000000];

s = 0;
% m_val is one of [3 5 7 9 13];
% theta_m for m=1:13.
theta = [%3.650024139523051e-008
  %5.317232856892575e-004
  1.495585217958292e-002  % m_vals = 3
  %8.536352760102745e-002
  2.539398330063230e-001  % m_vals = 5
  %5.414660951208968e-001
  9.504178996162932e-001  % m_vals = 7
  %1.473163964234804e+000
  2.097847961257068e+000  % m_vals = 9
  %2.811644121620263e+000
  %3.602330066265032e+000
  %4.458935413036850e+000
  5.371920351148152e+000];% m_vals = 13

T2 = T*T;
T4 = T2*T2;
T6 = T2*T4;
d4 = norm(T4,1)^(1/4);
d6 = norm(T6,1)^(1/6);
eta1 = max(d4, d6);
if (eta1 <= theta(1) && ell(T, coeff(1), 3) == 0)
  m = 3;
  return;
end
if (eta1 <= theta(2) && ell(T, coeff(2), 5) == 0)
  m = 5;
  return;
end

isSmall = size(T,1) < 150; %Compute matrix power explicitly
if isSmall
  d8 = norm(T4*T4,1)^(1/8);
else
  d8 = normAm(T4, 2)^(1/8);
end
eta3 = max(d6, d8);
if (eta3 <= theta(3) && ell(T, coeff(3), 7) == 0)
  m = 7;
  return;
end
if (eta3 <= theta(4) && ell(T, coeff(4), 9) == 0)
  m = 9;
  return;
end
if isSmall
  d10 = norm(T4*T6,1)^(1/10);
else
  d10 = normAm(T2, 5)^(1/10);
end
eta4 = max(d8, d10);
eta5 = min(eta3, eta4);
s = max(ceil(log2(eta5/theta(5))), 0);
s = s + ell(T/2^s, coeff(5), 13);
if isinf(s)
  % Overflow in ell subfunction. Revert to old estimate.
  [t, s] = log2(norm(T,1)/theta(end));
  s = s - (t == 0.5); % adjust s if normA/theta(end) is a power of 2.
end
m = 13;
end

function [c,mv] = normAm(A,m)
%NORMAM   Estimate of 1-norm of power of matrix.
%   NORMAM(A,M) estimates norm(A^m,1). If A has nonnegative elements
%   the estimate is exact.
%   [C, MV] = NORMAM(A,M) returns the estimate C and the number MV of
%   matrix-vector products computed involving A or A^*.

%   Reference:
%   A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring Algorithm
%      for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
%      970-989, 2009.
%
%   Awad H. Al-Mohy and Nicholas J. Higham
%   Copyright 2014-2015 The MathWorks, Inc.

n = size(A,1);
if n < 50 % Compute matrix power explicitly
  mv = 0;
  c = norm(matlab.internal.math.mpower.viaMtimes(A,m),1);
elseif isreal(A) && all(A(:) >= 0)
  % For positive matrices only.
  e = ones(n,1,class(A));
  for j=1:m
    e = A'*e;
  end
  c = norm(e,inf);
  mv = m;
else
  [c,~,~,it] = normest1(@afun_power);
  mv = it(2)*2*m; % Since t = 2.
end
% End of normAm

  function Z = afun_power(flag,X)
    %afun_power  Function to evaluate matrix products needed by normest1.
    if isequal(flag,'dim')
      Z = n;
    elseif isequal(flag,'real')
      Z = isreal(A);
    else
      if isequal(flag,'notransp')
        for i = 1:m
          X = A*X;
        end
      elseif isequal(flag,'transp')
        for i = 1:m
          X = A'*X;
        end
      end
      Z = X;
    end
  end
end

function D2 = wigD(q)
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
% D1 = zeros(3,3,size(q,2),size(q,3));
% D1(1,1,:,:) = A.^2;
% D1(1,2,:,:) = sqrt(2)*A.*B;
% D1(1,3,:,:) = B.^2;
% D1(2,1,:,:) = -sqrt(2)*A.*Bst;
% D1(2,2,:,:) = Z;
% D1(2,3,:,:) = sqrt(2)*Ast.*B;
% D1(3,1,:,:) = Bst.^2;
% D1(3,2,:,:) = -sqrt(2).*Ast.*Bst;
% D1(3,3,:,:) = Ast.^2;

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