% cardamom_propagatedm Propagate the density matrix of a spin-1/2 14N
%                      nitroxide using different methods from the literature.
%
%   Sprho = cardamom_propagatedm(Sys,Par,Opt,MD,omega,CenterField);
%
%     omega          double
%                    microwave frequency for CW field sweep, in Hz
%
%     CenterField    double
%                    magnetic field, in mT
%
%
%   Sys: stucture with system's dynamical parameters
%     g              numeric, size = (1,3)
%                    principal values of the g-tensor
%     A              numeric, size = (1,3)
%                    principal values of the A-tensor
%
%   Par: structure with simulation parameters
%     dt             double
%                    rotational dynamics propagation time step (in seconds)
%
%     Dt             double
%                    spin dynamics propagation time step (in seconds)
%
%   Exp: experimental parameter settings
%     B              double  TODO can this be replaced by a fieldsweep and then used to extract omega0?
%                    center magnetic field
%
%   Opt: optional settings
%     Method         string
%                    'Nitroxide': propagate using the m_s=-1/2 subspace
%                    'ISTOs': propagate using correlation functions
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

function Sprho = cardamom_propagatedm(Sys, Par, Opt, MD, omega, CenterField)

% Preprocessing
% -------------------------------------------------------------------------


persistent cacheTensors
persistent D2TrajMol
persistent gTensorState
persistent ATensorState
% persistent fullSteps
% persistent K2

if ~isfield(Opt, 'Liouville')
  if Opt.truncate
    % Liouville space propagation is needed for the ensemble-averaged
    % propagator
    Opt.Liouville = 1;
  else
    Opt.Liouville = 0;
  end
end

Liouville = Opt.Liouville;
Method = Opt.Method;
Model = Par.Model;

% define a rotational dynamics time scale for integrating corrlation
% functions
if ~isempty(MD)
  tauR = MD.tauR;
  isHMMfromMD = strcmp(MD.TrajUsage,'Markov');
  useMD = 1;
  isExplicit = strcmp(MD.TrajUsage,'Explicit');
else
  useMD = 0;
  isHMMfromMD = 0;
  if isfield(Sys,'Diff')       % TODO make this work for jumps and ISTOs
    tauR = 1/6/mean(Sys.Diff);
  end
end

if ~isfield(Par,'RTraj') && ~isfield(Par,'qTraj')
  error('Either Par.RTraj or Par.qTraj must be provided.')
end
% if ~isfield(MD,'RTraj') && (~isfield(Par,'RTraj')||~isfield(Par,'qTraj'))
%   error('Either Par.RTraj and Par.qTraj, or MD.RTraj must be provided.')
% end

% grab the quaternions or rotation matrices
if useMD
  if isExplicit
    if strcmp(Method, 'ISTOs')
      qTraj = Par.qTraj;
      qLab = Par.qLab;
    else
      RTraj = Par.RTraj;
      RTrajInv = permute(RTraj,[2,1,3,4]);
      RLab = Par.RLab;
      RLabInv = permute(RLab, [2,1,3,4]);
    end
  else
    if strcmp(Method, 'ISTOs')
      qTraj = Par.qTraj;
      qLab = Par.qLab;
    else
      RTraj = Par.RTraj;
      RTrajInv = permute(RTraj,[2,1,3,4]);
      RLab = Par.RLab;
      RLabInv = permute(RLab, [2,1,3,4]);
    end
  end
else
  if strcmp(Method, 'ISTOs')
    qTraj = Par.qTraj;
    qLab = Par.qLab;
  else
    RTraj = Par.RTraj;
    RTrajInv = permute(RTraj,[2,1,3,4]);
    RLab = Par.RLab;
    RLabInv = permute(RLab,[2,1,3,4]);
  end
end

if ~isfield(Sys,'g')
  error('g-tensor not specified.');
else
  g = Sys.g;
end

if isfield(Sys, 'A')
  A = Sys.A; 
else
  if strcmp(Method, 'Nitroxide')
    error('An A-tensor is required for the Nitroxide method.'); 
  end
end

Dt = Par.Dt;  % quantum spin propagation time step
nTraj = Par.nTraj;
isBlockAveraging = Par.isBlock;  % tensor time block averaging
nSteps = Par.nSteps;

truncate = Opt.truncate;  % ensemble averaged propagation

% store parameters for feeding into propagation function
Sim.nSteps = nSteps;
Sim.nTraj = nTraj;
Sim.truncate = truncate;
Sim.Liouville = Liouville;

% Gamma = gfree*bmagn/(planck/2/pi);  % rad s^-1 T^-1
% omegaN = 19.331e6*B;  % gyromagnetic ratio for 14N: 
                      % 19.331x10^6 rad s^-1 T^-1

% Simulation
% -------------------------------------------------------------------------

switch Method
  case 'Nitroxide'  % see Ref [1]
    
    if ~isequal(size(g),[1,3])
      error('g-tensor must be a 3-vector.')
    end

    if isfield(Sys, 'A')
      if ~isequal(size(A),[1,3])
        error('A-tensor must be a 3-vector.')
      end
    end
    
    if ~isHMMfromMD
      % Calculate time-dependent tensors from orientational trajectories
      gTensor = cardamom_tensortraj(g,RTraj,RTrajInv);
      ATensor = cardamom_tensortraj(A,RTraj,RTrajInv)*1e6*2*pi; % MHz (s^-1) -> Hz (rad s^-1)
      
    else
      % Calculate time-dependent tensors from state trajectories

      % Calculate the average interaction tensors for each state using 
      % MD-derived frame trajectories and Viterbi trajectories
      if isempty(gTensorState)
        % Perform MD-derived rotations on g- and A-tensors
        gTensorMD = cardamom_tensortraj(g,RTraj,RTrajInv);
        ATensorMD = cardamom_tensortraj(A,RTraj,RTrajInv)*1e6*2*pi; % MHz (s^-1) -> Hz (rad s^-1)
        
        nVitTraj = size(MD.viterbiTraj,1);
        gTensorState = zeros(3,3,MD.nStates,nVitTraj);
        ATensorState = zeros(3,3,MD.nStates,nVitTraj);

        % Average over time axis
        for iState = 1:MD.nStates
          for iTraj = 1:nVitTraj
            idxState = MD.viterbiTraj(iTraj,:) == iState;
            gTensorState(:,:,iState,iTraj) = squeeze(mean(gTensorMD(:,:,iTraj,idxState),4));
            ATensorState(:,:,iState,iTraj) = squeeze(mean(ATensorMD(:,:,iTraj,idxState),4));
          end
        end
        
        % Average over trajectories
        gTensorState = mean(gTensorState,4,'omitnan');
        ATensorState = mean(ATensorState,4,'omitnan');
      end
      
      % Calculate new time-dependent tensors from state trajectories
      % generated using optimized HMM parameters
      gTensor = zeros(3,3,nTraj,nSteps);
      ATensor = zeros(3,3,nTraj,nSteps);
      
      for iStep = 1:nSteps
        for iTraj = 1:nTraj
          state = Par.stateTraj(iTraj,iStep);
          gTensor(:,:,iTraj,iStep) = gTensorState(:,:,state);
          ATensor(:,:,iTraj,iStep) = ATensorState(:,:,state);
        end
      end
    end
    
    % Time block averaging and sliding window processing of tensors
    if isBlockAveraging
      gTensorBlock = zeros(3,3,size(gTensor,3),Par.nBlocks);
      ATensorBlock = zeros(3,3,size(ATensor,3),Par.nBlocks);

      % Average the interaction tensors over time blocks
      idx = 1:Par.BlockLength;
      for k = 1:Par.nBlocks
        gTensorBlock(:,:,:,k) = mean(gTensor(:,:,:,idx),4);
        ATensorBlock(:,:,:,k) = mean(ATensor(:,:,:,idx),4);
        idx = idx + Par.BlockLength;
      end

      if useMD && isExplicit
        % Perform sliding window processing if using MD trajectory explicitly
        gTensor = zeros(3,3,nTraj,nSteps);
        ATensor = zeros(3,3,nTraj,nSteps);
        for k = 1:nTraj
          idx = (1:nSteps) + (k-1)*Par.lag;
          gTensor(:,:,k,:) = gTensorBlock(:,:,idx);
          ATensor(:,:,k,:) = ATensorBlock(:,:,idx);
        end
      else
        % No sliding windows for other methods, as the length of generated
        % trajectories is determined by length of FID
        gTensor = gTensorBlock;
        ATensor = ATensorBlock;
      end
    end
    
    % Rotate tensors into lab frame explicitly
    if ~isempty(RLab)
      gTensor = multimatmult(RLab, multimatmult(gTensor, RLabInv));
      ATensor = multimatmult(RLab, multimatmult(ATensor, RLabInv));
    end
    
    gIso = sum(g)/3;
    GpTensor = (gTensor - gIso)/gfree;
    
    % Prepare propagators
    % ---------------------------------------------------------------------
    
    Gp_zz = GpTensor(3,3,:,:);

    % Norm of expression in Eq. 24 in [1]
    a = sqrt(ATensor(1,3,:,:).*ATensor(1,3,:,:) ...
           + ATensor(2,3,:,:).*ATensor(2,3,:,:) ...
           + ATensor(3,3,:,:).*ATensor(3,3,:,:));

    % Rotation angle and unit vector parallel to axis of rotation
    % refer to paragraph below Eq. 37 in [1]
    theta = Dt*0.5*squeeze(a);
    nx = squeeze(ATensor(1,3,:,:)./a);
    ny = squeeze(ATensor(2,3,:,:)./a);
    nz = squeeze(ATensor(3,3,:,:)./a);
    
    % Eqs. A1-A2 in [1]
    ct = cos(theta) - 1;
    st = -sin(theta);
    
    % Matrix exponential of hyperfine part
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
    U = bsxfun(@times, exp(-1i*Dt*0.5*omega*Gp_zz), expadotI);
    
    % Propagate density matrix
    % ---------------------------------------------------------------------
    
    rho = zeros(3,3,nTraj,nSteps);
    rho(:,:,:,1) = 0.5*repmat(eye(3),[1,1,nTraj]);

    Sim.fullSteps = nSteps;
    Sprho = propagate(rho, U, Method, Sim, Opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'ISTOs'  % see Ref [2]
    
    % Calculate and store rotational basis operators
    % ---------------------------------------------------------------------
    
    if isempty(cacheTensors)
      % ISTOs in the lab frame and IST components in the principal frame 
      % are time-independent, so we only need to calculate them once

      for iSpin = 1:numel(Sys.Spins)
        SpinOps{iSpin,1} = sop(Sys.Spins,iSpin,1);
        SpinOps{iSpin,2} = sop(Sys.Spins,iSpin,2);
        SpinOps{iSpin,3} = sop(Sys.Spins,iSpin,3);
      end

      % electron spin operators
      SpinOps{1,1} = zeros(size(SpinOps{1,1}));  % S_x and S_y are zero operators in HF limit TODO: generalize this
      SpinOps{1,2} = zeros(size(SpinOps{1,1}));

      if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end  % TODO include frames in cardamom

      [T,F,~,~,~] = magint(Sys,SpinOps,CenterField,0,0);

%       F0 = F.F0*2*pi;
      if isfield(Sys, 'A')
        F0 = F.F0(2)*2*pi;  % Hz -> rad s^-1, only keep isotropic HF interaction
      end
      F2 = F.F2*2*pi;  % Hz -> rad s^-1

%       T0 = T.T0;
      if isfield(Sys, 'A')
        T0 = T.T0{2};  % only keep isotropic HF interaction
      end
      T2 = T.T2;

      % zeroth rank
  %       cacheTensors.Q0 = conj(F0(1))*T0{1} + conj(F0(2))*T0{2};

      if isfield(Sys, 'A')
%         if Liouville
%           cacheTensors.Q0 = conj(F0)*tosuper(T0,'c');
%         else
          cacheTensors.Q0 = conj(F0)*T0;
%         end
      else
        cacheTensors.Q0 = 0;
      end

      cacheTensors.Q2 = cell(5,5);

      % create the 25 second-rank RBOs
      for mp = 1:5
        for m = 1:5
%           if Liouville
%             cacheTensors.Q2{mp,m} = zeros(size(T2{1}).^2);
%             for iSpin = 1:numel(Sys.Spins)
%               cacheTensors.Q2{mp,m} = cacheTensors.Q2{mp,m} ...
%                                       + conj(F2(iSpin,mp))*tosuper(T2{iSpin,m},'c');
%             end
%           else
            cacheTensors.Q2{mp,m} = zeros(size(T2{1}));
            for iSpin = 1:numel(Sys.Spins)
              cacheTensors.Q2{mp,m} = cacheTensors.Q2{mp,m} + conj(F2(iSpin,mp))*T2{iSpin,m};
            end
%           end
        end
      end

    end
    
    % Process Wigner D-matrices
    % ---------------------------------------------------------------------

    % time block averaging and sliding window processing
    if isBlockAveraging
      if isempty(D2TrajMol)
        % if using an MD trajectory, it is best to process the molecular
        % dynamics once and store the result
        % this variable will be set to empty later if an MD trajectory is
        % not being used
        D2TrajMol = wigD(qTraj);
        D2Avg = zeros(5,5,Par.nBlocks);

        % average the Wigner D-matrices over time blocks
        idx = 1:Par.BlockLength;
        for k = 1:Par.nBlocks
          D2Avg(:,:,k) = mean(D2TrajMol(:,:,:,idx),4);
          idx = idx + Par.BlockLength;
        end

        % perform sliding window processing if using MD trajectory explicitly
        if useMD
          if isExplicit
            D2TrajMol = zeros(5,5,nTraj,nSteps);

            for k = 1:nTraj
              idx = (1:nSteps) + (k-1)*Par.lag;
              D2TrajMol(:,:,k,:) = D2Avg(:,:,idx);
            end
          else
            D2TrajMol = D2Avg;
          end
        end
      end
      D2Traj = D2TrajMol;
      if ~useMD
        % trajectories will continually be generated, so leave this empty
        D2TrajMol = [];
      end
    else
      D2Traj = wigD(qTraj);
    end
    
    % check for new lab frame rotations
    if ~isempty(qLab)
      D2Lab = wigD(qLab);
      D2Traj = multimatmult(D2Lab, D2Traj);
    end
    
    if truncate
%       if isempty(fullSteps)||isempty(K2)

      % set the number of time steps to use full propagator before 
      % correlation functions relax
      if strcmp(truncate, 'all')
        fullSteps = 1;  % use ensemble-averaged propagator the entire time
      else
        fullSteps = ceil(truncate*1e-9/Dt);
      end
      
      % if correlation time is longer than user input, just use full 
      % propagation scheme the entire time
      if fullSteps>nSteps
        fullSteps = nSteps; 
        truncate = 0;
      else

        % approximate upper limit of correlation function integration based
        % on rotational correlation time
        NtauR = 10*ceil(tauR/Dt);
        
        if NtauR>nSteps
          NtauR = nSteps;
        end
        
        % calculate correlation functions for approximate propagator
        D2Acorr = autocorrfft(D2Traj, 4, 0, 0, 1);
%         D2Acorr = autocorrfft(D2Traj-mean(D2Traj,4), 4, 0, 0, 0);
  %       for mppp=1:5
  %       for mpp=1:5
  %       for mp=1:5
  %         for m=1:5
  %           % non-normalized autocorrelation functions are needed here
  %           D2Traj1 = squeeze(D2(mp,m,:,:));
  %           D2Traj2 = squeeze(D2(mppp,mpp,:,:));
  %           xcorr = crosscorrfft(D2Traj2, D2Traj1, 2, 0, 0, 0);
  %           xcorrAvg = mean(xcorr,1);
  %           K2(mppp,mpp,mp,m) = trapz(time(1:NtauR), xcorrAvg(1:NtauR));
  %         end
  %       end
  %       end
  %       end

        if strcmp(Opt.debug.EqProp,'time')
          K2 = squeeze(trapz(D2Acorr(:,:,:,1:NtauR),4))*Dt;  % TODO implement cross-correlation functions
        elseif strcmp(Opt.debug.EqProp,'all')
          K2 = squeeze(trapz(mean(D2Acorr(:,:,:,1:NtauR),3),4))*Dt;  % TODO implement cross-correlation functions
        end

        idx = K2>1e-11;
        K2 = K2.*idx;

      end
%       end

    else
      % no approximation, use full propagator for all time steps
      fullSteps = nSteps;

    end

    D2avg = mean(D2Traj, 4);  % time average of trajectories of Wigner matrices
    
    % calculate Hamiltonians and propagators
    % ---------------------------------------------------------------------
    
    H = repmat(cacheTensors.Q0,[1,1,nTraj,fullSteps]);

    % rotate second rank terms and add to Hamiltonian
    for mp = 1:5
      for m = 1:5
        H = H + bsxfun(@times, D2Traj(m,mp,:,1:fullSteps), cacheTensors.Q2{mp,m});
      end
    end
    
%     if Liouville
%       H = tosuper(H, 'c');
%     end
    
    U = zeros(size(H));
    
    for iStep=1:fullSteps
      for iTraj=1:nTraj
        U(:,:,iTraj,iStep) = expeig(1i*Dt*H(:,:,iTraj,iStep));  % TODO speed this up!
%         U(:,:,iTraj,iStep) = expm_fast1(1i*dt*H(:,:,iTraj,iStep));
      end
    end

    if truncate
      % Prepare equilibrium propagators
      % -------------------------------------------------------------------

%       if Liouville
%         Heq = repmat(cacheTensors.Q0,[1,1,nTraj]);
%       else
        Heq = repmat(tosuper(cacheTensors.Q0,'c'),[1,1,nTraj]);
%       end

      HeqOrder1 = 0;
      HeqOrder2 = 0;

      % rotate second rank terms and add to Hamiltonian
      for mp = 1:5
        for m = 1:5
%           if Liouville
%             HeqOrder1 = HeqOrder1 + bsxfun(@times, D2avg(m,mp,:), cacheTensors.Q2{mp,m});
%             HeqOrder2 = HeqOrder2 + 1i*cacheTensors.Q2{mp,m}*(cacheTensors.Q2{mp,m})'.*K2(mp,m,:);
%           else
            HeqOrder1 = HeqOrder1 + bsxfun(@times, D2avg(m,mp,:), tosuper(cacheTensors.Q2{mp,m},'c'));
            HeqOrder2 = HeqOrder2 + 1i*tosuper(cacheTensors.Q2{mp,m},'c')*tosuper(cacheTensors.Q2{mp,m}','c').*K2(mp,m,:);
%           end

%       for mppp = 1:5
    %         for mpp = 1:5
    %           HeqOrder2 = HeqOrder2 + 1i*cacheTensors.Q2{mp,m}*(cacheTensors.Q2{mppp,mpp})'*K2(mppp,mpp,mp,m,:);
    %         end
    %       end
        end
      end

      Heq = bsxfun(@plus,Heq+HeqOrder1,HeqOrder2);  % FIXME plus or minus?
%       Heq = Heq + HeqOrder2;  % FIXME plus or minus?

%       if Liouville
        if strcmp(Opt.debug.EqProp,'time')
          Ueq = zeros(size(U,1)^2,size(U,2)^2,nTraj);
          for iTraj=1:nTraj
            Ueq(:,:,iTraj) = expm(1i*Dt*Heq(:,:,iTraj));
          end
        elseif strcmp(Opt.debug.EqProp,'all')
          Ueq = expm(1i*Dt*mean(Heq,3));
        end
        Sim.Ueq = Ueq;
%       else
%         Ueq = zeros(size(U,1),size(U,2),nTraj);
%         Ueqdag = zeros(size(U,1),size(U,2),nTraj);
%         for iTraj=1:nTraj
%           Ueq(:,:,iTraj) = expm_fast1(1i*dt*Heq(:,:,iTraj));
%           Ueqdag(:,:,iTraj) = expm_fast1(-1i*dt*Heq(:,:,iTraj));
%         end
%         Sim.Ueq = Ueq;
%         Sim.Ueqdag = Ueqdag;
%       end

%       Ueqdag = conj(permute(Ueq, [2,1,3,4]));

    end

    % Set up starting state of density matrix after pi/2 pulse
    % ---------------------------------------------------------------------

    Sim.fullSteps = fullSteps;
    rho = zeros(Sys.nStates,Sys.nStates,nTraj,nSteps);
    
    % initial state after pi/2 pulse is |S_x>|E_I>
    if isfield(Sys, 'A')
      % 1 electron + 1 nucleus
      rho0 = sop(Sys.Spins,'xe');
    else
      % 1 electron only
      rho0 = sop(Sys.S,'x');
    end
    
    rho(:,:,:,1) = repmat(rho0,1,1,nTraj,1);
    
    % Propagate density matrix
    % ---------------------------------------------------------------------
    
    rho = propagate(rho, U, Method, Sim, Opt);
    
    if Liouville
      if strcmp(Opt.debug.EqProp,'time')
        rho = reshape(rho,[Sys.nStates,Sys.nStates,nTraj,nSteps]);
      elseif strcmp(Opt.debug.EqProp,'all')
        rho = reshape(rho,[Sys.nStates,Sys.nStates,nSteps]);
      end
    end
    
    % Multiply density matrix result by S_+ detection operator
    % ---------------------------------------------------------------------
    
    % Only keep the m_S=\beta subspace part that contributes to 
    %   tr(S_{+}\rho(t))
    if isfield(Sys, 'A')
      % 1 electron + 1 nucleus
      projector = sop(Sys.Spins,'+e');
    else
      % 1 electron only
      projector = sop(Sys.S,'+');
    end
    
    if strcmp(Opt.debug.EqProp,'time')
      Sprho = multimatmult(repmat(projector,1,1,nTraj,nSteps),rho);
    elseif strcmp(Opt.debug.EqProp,'all')
      Sprho = multimatmult(repmat(projector,1,1,nSteps),rho);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise
    error('Propagation method not recognized.')
end


end

% Helper functions
% -------------------------------------------------------------------------

function rho = propagate(rho, U, Method, Sim, Opt)

fullSteps = Sim.fullSteps;
nSteps = Sim.nSteps;
nTraj = Sim.nTraj;
truncate = Sim.truncate;
Liouville = Sim.Liouville;

% Propagate density matrix
% ---------------------------------------------------------------------

% N = size(U,1);

% if Liouville
%   rho = reshape(rho,[N,1,nTraj,nSteps]);
% end

% propagation using full Hamiltonian
switch Method
  case 'Nitroxide'
    % subspace propagation only requires U, but not Udag
    for iStep=2:fullSteps
      rho(:,:,:,iStep) = multimatmult( U(:,:,:,iStep-1),...
                              multimatmult( rho(:,:,:,iStep-1),...
                                   U(:,:,:,iStep-1) ) );                  
    end

  case 'ISTOs'
    
%     if Liouville
%       for iStep=2:fullSteps
%         rho(:,:,:,iStep) = multimatmult(U(:,:,:,iStep-1), rho(:,:,:,iStep-1));   
% %         for iTraj=1:nTraj
% %           rho(:,:,iTraj,iStep) = U(:,:,iTraj,iStep-1)*rho(:,:,iTraj,iStep-1);  
% %         end
%       end

%     else % Hilbert space
      Udag = conj(permute(U,[2,1,3,4]));
      for iStep=2:fullSteps
        rho(:,:,:,iStep) = multimatmult( U(:,:,:,iStep-1),...
                                   multimatmult( rho(:,:,:,iStep-1),...
                                         Udag(:,:,:,iStep-1) ) );                  
      end
%     end
end

if strcmp(Opt.debug.EqProp,'all')
  rho = squeeze(mean(rho, 3));
end

if truncate
  % propagation using correlation function approximation
%   if Liouville
    Ueq = Sim.Ueq;
    if strcmp(Opt.debug.EqProp,'time')
      rho = reshape(rho,[size(Ueq,1),1,nTraj,nSteps]);
      for iStep=fullSteps+1:nSteps
        rho(:,:,:,iStep) = multimatmult(Ueq, rho(:,:,:,iStep-1)); 
      end
    elseif strcmp(Opt.debug.EqProp,'all')
      rho = reshape(rho,[size(Ueq,1),nSteps]);
      for iStep=fullSteps+1:nSteps
        rho(:,iStep) = Ueq*rho(:,iStep-1); 
      end
    end

%   else
%     Ueq = Sim.Ueq;
%     Ueqdag = Sim.Ueqdag;
%     for iStep=fullSteps+1:nSteps
%       rho(:,:,:,iStep) = multimatmult( Ueq, multimatmult( rho(:,:,:,iStep-1), Ueqdag ) );  
%     end
%   end

end

end

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

ASq = A.*A;
BSq = B.*B;
AstSq = Ast.*Ast;
BstSq = Bst.*Bst;

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
% D2(1,1,:,:) = A.^4;
D2(1,1,:,:) = ASq.*ASq;
D2(1,2,:,:) = 2*A.^3.*B;
% D2(1,3,:,:) = sqrt(6)*A.^2.*B.^2;
D2(1,3,:,:) = sqrt(6)*ASq.*BSq;
D2(1,4,:,:) = 2*A.*B.^3;
% D2(1,5,:,:) = B.^4;
D2(1,5,:,:) = BSq.*BSq;
D2(2,1,:,:) = -2*A.^3.*Bst;
% D2(2,2,:,:) = A.^2.*(2*Z-1);
D2(2,2,:,:) = ASq.*(2*Z-1);
D2(2,3,:,:) = sqrt(6)*A.*B.*Z;
% D2(2,4,:,:) = B.^2.*(2*Z+1);
D2(2,4,:,:) = BSq.*(2*Z+1);
D2(2,5,:,:) = 2*Ast.*B.^3;
% D2(3,1,:,:) = sqrt(6)*A.^2.*Bst.^2;
D2(3,1,:,:) = sqrt(6)*ASq.*BstSq;
D2(3,2,:,:) = -sqrt(6)*A.*Bst.*Z;
D2(3,3,:,:) = 1/2*(3*Z.^2-1);
D2(3,4,:,:) = sqrt(6)*Ast.*B.*Z;
% D2(3,5,:,:) = sqrt(6)*Ast.^2.*B.^2;
D2(3,5,:,:) = sqrt(6)*AstSq.*BSq;
D2(4,1,:,:) = -2*A.*Bst.^3;
% D2(4,2,:,:) = Bst.^2.*(2*Z+1);
D2(4,2,:,:) = BstSq.*(2*Z+1);
D2(4,3,:,:) = -sqrt(6)*Ast.*Bst.*Z;
% D2(4,4,:,:) = Ast.^2.*(2*Z-1);
D2(4,4,:,:) = AstSq.*(2*Z-1);
D2(4,5,:,:) = 2*Ast.^3.*B;
% D2(5,1,:,:) = Bst.^4;
D2(5,1,:,:) = BstSq.*BstSq;
D2(5,2,:,:) = -2*Ast.*Bst.^3;
% D2(5,3,:,:) = sqrt(6)*Ast.^2.*Bst.^2;
D2(5,3,:,:) = sqrt(6)*AstSq.*BstSq;
D2(5,4,:,:) = -2*Ast.^3.*Bst;
% D2(5,5,:,:) = Ast.^4;
D2(5,5,:,:) = AstSq.*AstSq;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deprecated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   case 'Resampling'
%   
%     g_tr = sum(g);
%     
%     % Process MD simulation data
%     % ---------------------------------------------------------------------
%     
%     % Perform rotations on g- and A-tensors
%     gTensor = tensortraj(g,RTraj,RTrajInv);
% 
%     ATensor = tensortraj(A,RTraj,RTrajInv)*1e6*2*pi; % MHz (s^-1) -> Hz (rad s^-1)
%     
%     if isfield(MD,'GlobalDiff')
%       gTensor = multimatmult(RTrajGlobal, multimatmult(gTensor, RTrajGlobalInv, 'real'), 'real');
%       ATensor = multimatmult(RTrajGlobal, multimatmult(ATensor, RTrajGlobalInv, 'real'), 'real');
%     end
%     
%     GpTensor = gTensor/gfree - g_tr/3/gfree;
%         
%     rho = zeros(3,3,Par.nTraj,nSteps);
%     rho(:,:,:,1) = repmat(eye(3),[1,1,Par.nTraj]);
%     
%     % Prepare propagators
%     % ---------------------------------------------------------------------
%     
%     Gp_zz = GpTensor(3,3,:,:);
% 
%     % norm of expression in Eq. 24 in [1]
%     a = sqrt(ATensor(1,3,:,:).*ATensor(1,3,:,:) ...
%            + ATensor(2,3,:,:).*ATensor(2,3,:,:) ...
%            + ATensor(3,3,:,:).*ATensor(3,3,:,:));
% 
%     % rotation angle and unit vector parallel to axis of rotation
%     % refer to paragraph below Eq. 37 in [1]
% %     theta = Gamma*dt*0.5*squeeze(a);
%     theta = dt*0.5*squeeze(a);
%     nx = squeeze(ATensor(1,3,:,:)./a);
%     ny = squeeze(ATensor(2,3,:,:)./a);
%     nz = squeeze(ATensor(3,3,:,:)./a);
% 
%     % Eqs. A1-A2 in [1]
%     ct = cos(theta) - 1;
%     st = -sin(theta);
% 
%     % matrix exponential of hyperfine part
%     % Eqs. A1-A2 are used to construct Eq. 37 in [1]
%     expadotI = zeros(3,3,size(Gp_zz,3),size(Gp_zz,4));
%     expadotI(1,1,:,:) = 1 + ct.*(nz.*nz + 0.5*(nx.*nx + ny.*ny)) ...
%                          + 1i*st.*nz;
%     expadotI(1,2,:,:) = sqrt(0.5)*(st.*ny + ct.*nz.*nx) ...
%                          + 1i*sqrt(0.5)*(st.*nx - ct.*nz.*ny);
%     expadotI(1,3,:,:) = 0.5*ct.*(nx.*nx - ny.*ny) ... 
%                          - 1i*ct.*nx.*ny;
%     expadotI(2,1,:,:) = sqrt(0.5)*(-st.*ny + ct.*nz.*nx) ...
%                          + 1i*sqrt(0.5)*(st.*nx + ct.*nz.*ny);
%     expadotI(2,2,:,:) = 1 + ct.*(nx.*nx + ny.*ny);
%     expadotI(2,3,:,:) = sqrt(0.5)*(st.*ny - ct.*nz.*nx) ...
%                          + 1i*sqrt(0.5)*(st.*nx + ct.*nz.*ny);
%     expadotI(3,1,:,:) = 0.5*ct.*(nx.*nx - ny.*ny) ...
%                          + 1i*ct.*nx.*ny;
%     expadotI(3,2,:,:) = sqrt(0.5)*(-st.*ny - ct.*nz.*nx) ...
%                          + 1i*sqrt(0.5)*(st.*nx - ct.*nz.*ny);
%     expadotI(3,3,:,:) = 1 + ct.*(nz.*nz + 0.5*(nx.*nx + ny.*ny))  ...
%                          - 1i*st.*nz;
%     
%     
%     % Calculate propagator, Eq. 35 in [1]
%     U = bsxfun(@times, exp(-1i*dt*0.5*omega*Gp_zz), expadotI);
%     
%     
%     % Propagate density matrix
%     % ---------------------------------------------------------------------
%     
%     if nTraj>1
%       for iStep=2:nSteps
%         rho(:,:,:,iStep) = multimatmult(U(:,:,:,iStep-1),...
%                                    multimatmult(rho(:,:,:,iStep-1), U(:,:,:,iStep-1),'complex'),...
%                                    'complex');
% 
%         %  Trace of a product of matrices is the sum of entry-wise products
%         %      rho_t(:,:,:,iStep) = U2(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
%       end
%     else
%       for iStep=2:nSteps
%         rho(:,:,1,iStep) = U(:,:,1,iStep-1)*rho(:,:,1,iStep-1)*U(:,:,1,iStep-1);
% 
%         %  Trace of a product of matrices is the sum of entry-wise products
%         %      rho_t(:,:,:,iStep) = U2(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
%       end
%     end

% else
% % truncated trajectory behavior not being used, so use full propagation
% % scheme for all time steps
%   
%   switch Method
%     case 'Nitroxide'
%       for iStep=2:nSteps
%         rho(:,:,:,iStep) = multimatmult(U(:,:,:,iStep-1),...
%                                    multimatmult(rho(:,:,:,iStep-1), U(:,:,:,iStep-1),'complex'),...
%                                    'complex');
%     
%         %  Trace of a product of matrices is the sum of entry-wise products
%         %      rho_t(:,:,:,iStep) = U2(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
%       end
%       
%     case 'ISTOs'
%       if Liouville
%         for iStep=2:nSteps
%           for iTraj=1:nTraj
%             rho(:,iTraj,iStep) = U(:,:,iTraj,iStep-1)*rho(:,iTraj,iStep-1);  
%           end
%         end
% 
%       else % Hilbert space
%         for iStep=2:nSteps
%           rho(:,:,:,iStep) = multimatmult(U(:,:,:,iStep-1),...
%                                      multimatmult(rho(:,:,:,iStep-1),...
%                                            Udag(:,:,:,iStep-1),'complex'),...
%                                      'complex');                  
%         end
%       end
%   end
% 
% end

% Density Matrix Propagation
%     if nTraj>1
%       for iStep=2:nSteps
%         rho_t(:,:,:,iStep) = multimatmult(U(:,:,:,iStep-1),...
%                                    multimatmult(rho_t(:,:,:,iStep-1), U(:,:,:,iStep-1),'complex'),...
%                                    'complex');
%     
%         %  Trace of a product of matrices is the sum of entry-wise products
%         %      rho_t(:,:,:,iStep) = U2(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
%       end
%     else
%       for iStep=2:nSteps
%         rho_t(:,:,1,iStep) = U(:,:,1,iStep-1)*rho_t(:,:,1,iStep-1)*U(:,:,1,iStep-1);
% 
%         %  Trace of a product of matrices is the sum of entry-wise products
%         %      rho_t(:,:,:,iStep) = U2(:,:,:,iStep-1).*rho_t(:,:,:,iStep-1);
%       end
%     end