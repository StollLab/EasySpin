% cardamom_propagatedm Propagate the density matrix of a spin system.
%
%   Sprho = cardamom_propagatedm(Sys,Par,Opt,MD,omega,CenterField);
%
% Inputs:
%   Sys: stucture with system's dynamical parameters
%     g              3-array, principal values of the g-tensor
%     A              3-array, principal values of the A-tensor
%   Par: structure with simulation parameters
%     dtSpin         spin dynamics propagation time step (in seconds)
%     nSteps         number of steps
%     BlockLength    block length for block averaging, no averaging if set to 1
%     Model          'MD-HMM','MD-direct','MD-HBD','diffusion','jump'
%     RTraj          rotation matrix trajectory
%   Exp: experimental parameter settings
%     B              center magnetic field
%   Opt: optional settings
%     Method         'fast': propagate using the m_s=-1/2 subspace
%                    'ISTOs': propagate using correlation functions
%   MD:
%     RTraj          numeric, size = (3,3,nTraj,nSteps)
%                    externally provided rotation matrices
%   omega            double
%                    microwave frequency for CW field sweep, in Hz
%   CenterField      double
%                    magnetic field, in mT
%
%  Output:
%    Sprho          numeric, size = (3,3,nTraj,nSteps)
%                   a series of density matrices

% Implementations based on 
% [1] Sezer, et al., J. Chem. Phys. 128, 165106 (2008)
%      http://dx.doi.org/10.1063/1.2908075 
% [2] Oganesyan, Phys. Chem. Chem. Phys. 13, 4724 (2011)
%      http://dx.doi.org/10.1039/c0cp01068e

function Sprho = cardamom_propagatedm(Sys, Par, Opt, MD, omega, CenterField)

persistent cacheTensors cacheD2Traj gTensorState ATensorState

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0;
end

global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

% Preprocessing
%-------------------------------------------------------------------------------

PropagationMethod = Opt.Method;
switch PropagationMethod
  case 'fast'
    if ~isfield(Par,'RTraj')
      error('Par.RTraj must be provided.');
    end
  case 'ISTOs'
    if ~isfield(Par,'qTraj')
      error('Par.qTraj must be provided.');
    end
  otherwise
    error('Unknown propagation method ''%s''.',PropagationMethod);
end

% Define a rotational dynamics time scale for integrating correlation functions
useMD = ~isempty(MD);
if useMD
  tauR = MD.tauR;
  isHMMfromMD = strcmp(Par.Model,'MD-HMM');
  isMDdirect = strcmp(Par.Model,'MD-direct');
else
  isHMMfromMD = false;
  isMDdirect = false;
  if isfield(Sys,'Diff')       % TODO make this work for jumps and ISTOs
    tauR = 1/6/mean(Sys.Diff);
  end
end


dtSpin = Par.dtSpin;  % spin dynamics propagation time step
nSteps = Par.nSteps; % number of spin dynamics propagation steps
nTraj = Par.nTraj;

doBlockAveraging = Par.BlockLength>1;  % tensor time block averaging
doSlidingWindowProcessing = isMDdirect && Par.lag>0; % sliding window processing
alternatingTrajDirection = false;

% Simulation
%-------------------------------------------------------------------------------

switch PropagationMethod
  case 'fast'  % see Ref [1]
    
    if ~isfield(Sys,'g')
      error('A g-tensor is required for the fast method.');
    end
    g = Sys.g;
    if ~isequal(size(g),[1,3])
      error('g-tensor must be a 3-vector.')
    end
    
    includeHF = isfield(Sys,'A');
    if includeHF
      A = Sys.A;
      if ~isequal(size(A),[1,3])
        error('A-tensor must be a 3-vector.')
      end
    end
        
    % Calculate tensor trajectories
    %---------------------------------------------------------------------------
    RTrajInv = permute(Par.RTraj,[2,1,3,4]);
    if ~isHMMfromMD
      logmsg(2,'  calculating tensor trajectories from orientational trajectories');
      % Calculate time-dependent tensors from orientational trajectories
      gTensor = cardamom_tensortraj(g,Par.RTraj,RTrajInv);
      if includeHF
        ATensor = cardamom_tensortraj(A,Par.RTraj,RTrajInv)*1e6*2*pi; % MHz (s^-1) -> Hz (rad s^-1)
      end
      
    else
      % Calculate time-dependent tensors from state trajectories
      
      % Calculate the average interaction tensors for each state using
      % MD-derived frame trajectories and Viterbi state trajectories
      if isempty(gTensorState)
        logmsg(2,'  calculating tensor trajectories from orientational trajectories');
        % Perform MD-derived rotations on g- and A-tensors
        gTensorMD = cardamom_tensortraj(g,Par.RTraj,RTrajInv);
        if includeHF
          ATensorMD = cardamom_tensortraj(A,Par.RTraj,RTrajInv)*1e6*2*pi; % MHz (s^-1) -> Hz (rad s^-1)
        end
        
        % Average over time axis
        logmsg(1,'  calculating state averages of tensors');
        nVitTraj = size(MD.viterbiTraj,1);
        gTensorState = zeros(3,3,MD.nStates,nVitTraj);
        if includeHF
          ATensorState = zeros(3,3,MD.nStates,nVitTraj);
        end
        for iState = 1:MD.nStates
          for iTraj = 1:nVitTraj
            idxState = MD.viterbiTraj(iTraj,:) == iState;
            gTensorState(:,:,iState,iTraj) = squeeze(mean(gTensorMD(:,:,idxState,iTraj),3));
            if includeHF
              ATensorState(:,:,iState,iTraj) = squeeze(mean(ATensorMD(:,:,idxState,iTraj),3));
            end
          end
        end
        
        % Average over trajectories
        gTensorState = mean(gTensorState,4,'omitnan');
        if includeHF
          ATensorState = mean(ATensorState,4,'omitnan');
        end
      end
      
      % Calculate new time-dependent tensors from state trajectories
      % generated using optimized HMM parameters
      logmsg(2,'  calculating tensor trajectories from state trajectories');
      gTensor = zeros(3,3,nSteps,nTraj);
      if includeHF
        ATensor = zeros(3,3,nSteps,nTraj);
      end      
      for iState = 1:MD.nStates
        for iTraj = 1:nTraj
          idxState = Par.stateTraj(:,iTraj)==iState;
          gTensor(:,:,idxState,iTraj) = repmat(gTensorState(:,:,iState),[1,1,sum(idxState),1]);
          if includeHF
            ATensor(:,:,idxState,iTraj) = repmat(ATensorState(:,:,iState),[1,1,sum(idxState),1]);
          end
        end
      end
      
    end
    
    % Coarse-grain tensor trajectories by block averaging
    %---------------------------------------------------------------------------
    if doBlockAveraging
      logmsg(2,'  coarse-graining with tensor averaging, block length %d frames',Par.BlockLength);
      nBlocks = floor(size(gTensor,3)/Par.BlockLength);
      gTensorBlock = zeros(3,3,nBlocks,size(gTensor,4));
      if includeHF
        ATensorBlock = zeros(3,3,nBlocks,size(ATensor,4));
      end
      idxBlock = 1:Par.BlockLength;
      for iBlock = 1:nBlocks
        gTensorBlock(:,:,iBlock,:) = mean(gTensor(:,:,idxBlock,:),3);
        if includeHF
          ATensorBlock(:,:,iBlock,:) = mean(ATensor(:,:,idxBlock,:),3);
        end
        idxBlock = idxBlock + Par.BlockLength;
      end
    else
      logmsg(2,'  no coarse-graining tensor averaging');
      gTensorBlock = gTensor;
      if includeHF
        ATensorBlock = ATensor;
      end
    end
    
    % Sliding window processing: generate multiple trajectories from single one
    %---------------------------------------------------------------------------
    if doSlidingWindowProcessing
      logmsg(2,'  sliding window processing, increment %d',Par.lag);
      gTensor = zeros(3,3,nSteps,nTraj);
      if includeHF
        ATensor = zeros(3,3,nSteps,nTraj);
      end
      idx = 1:nSteps;
      for iTraj = 1:nTraj
        gTensor(:,:,:,iTraj) = gTensorBlock(:,:,idx,1);
        if includeHF
          ATensor(:,:,:,iTraj) = ATensorBlock(:,:,idx,1);
        end
        idx = idx + Par.lag;
        if alternatingTrajDirection
          idx = flip(idx);
        end
      end
    else
      logmsg(2,'  no sliding window processing');
      gTensor = gTensorBlock;
      if includeHF
        ATensor = ATensorBlock;
      end
    end
    
    % Combine local and global dynamics
    %---------------------------------------------------------------------------
    logmsg(2,'  combine local with global dynamics');
    if ~isempty(Par.RLab)
      RLabInv = permute(Par.RLab,[2,1,3,4]);
      gTensor = multimatmult(Par.RLab, multimatmult(gTensor, RLabInv));
      if includeHF
        ATensor = multimatmult(Par.RLab, multimatmult(ATensor, RLabInv));
      end
    end
    
    % Prepare propagators
    %---------------------------------------------------------------------------
    gIso = sum(g)/3;
    GpTensor = (gTensor - gIso)/gfree;    
    Gp_zz = GpTensor(3,3,:,:);

    if includeHF
      % Norm of expression in Eq. 24 in [1]
      a = sqrt(ATensor(1,3,:,:).*ATensor(1,3,:,:) ...
        + ATensor(2,3,:,:).*ATensor(2,3,:,:) ...
        + ATensor(3,3,:,:).*ATensor(3,3,:,:));
      
      % Rotation angle and unit vector parallel to axis of rotation
      % refer to paragraph below Eq. 37 in [1]
      theta = dtSpin*0.5*squeeze(a);
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
      U = bsxfun(@times, exp(-1i*dtSpin*0.5*omega*Gp_zz), expadotI);
    else
      U = exp(-1i*dtSpin*0.5*omega*Gp_zz);
    end
    
    % Set up starting state of density matrix after pi/2 pulse, S_x
    %---------------------------------------------------------------------------
    if includeHF
      rho = zeros(3,3,nSteps,nTraj);
      rho(:,:,1,:) = 0.5*repmat(eye(3),[1,1,nTraj]);
    else
      rho = zeros(1,1,nSteps,nTraj);
      rho(:,:,1,:) = 0.5;
    end
    
    % Propagate density matrix
    %---------------------------------------------------------------------------
    logmsg(2,'  propagate density matrix');
    Sprho = propagate(rho,U,nSteps,PropagationMethod);
    
    % Average over trajectories (4th dimension)
    Sprho = mean(Sprho,4);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'ISTOs'  % see Ref [2]
    
    % Calculate and cache rotational basis operators
    %---------------------------------------------------------------------------
    if isempty(cacheTensors)
      % ISTOs in the lab frame and IST components in the principal frame 
      % are time-independent, so we only need to calculate them once
      
      for iSpin = 1:numel(Sys.Spins)
        SpinOps{iSpin,1} = sop(Sys.Spins,[iSpin,1]);
        SpinOps{iSpin,2} = sop(Sys.Spins,[iSpin,2]);
        SpinOps{iSpin,3} = sop(Sys.Spins,[iSpin,3]);
      end
      
      if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end  % TODO include frames in cardamom
      
      [T,F,~,~,~] = magint(Sys,SpinOps,CenterField,false,false,false);
      F0 = F.F0*2*pi;  % Hz -> rad s^-1
      F2 = F.F2*2*pi;
      T0 = T.T0;
      T2 = T.T2;
      
      % zeroth rank
      cacheTensors.Q0 = 0;
      for k = 1:numel(F0)
        cacheTensors.Q0 = cacheTensors.Q0 + conj(F0(k))*T0{k};
      end
      
      % isotropic Zeeman interaction, for interaction frame transformation
      cacheTensors.H0 = conj(F0(1))*T0{1};
      
      % second rank (25 "RBOs")
      cacheTensors.Q2 = cell(5,5);
      for mp = 1:5
        for m = 1:5
          cacheTensors.Q2{mp,m} = zeros(size(T2{1}));
          for iInt = 1:size(F2,1)
            cacheTensors.Q2{mp,m} = cacheTensors.Q2{mp,m} + conj(F2(iInt,mp))*T2{iInt,m};
          end
        end
      end
      
    end
    
    % D2 trajectories, incl. time block averaging and sliding window processing
    %---------------------------------------------------------------------------
    if ~isempty(cacheD2Traj) % use cached D2 trajectories if present
      D2Traj = cacheD2Traj;
    else % if not cached, process q trajectory to get D2 trajectories
      
      D2Traj = wigD(Par.qTraj);

      % Average the Wigner D-matrices over time blocks
      if doBlockAveraging
        nBlocks = floor(size(D2Traj,3)/Par.BlockLength);
        D2BlockAvg = zeros(5,5,nBlocks,size(D2Traj,4));
        idx = 1:Par.BlockLength;
        for iBlock = 1:nBlocks
          D2BlockAvg(:,:,iBlock,:) = mean(D2Traj(:,:,idx,:),3);
          idx = idx + Par.BlockLength;
        end
      else
        D2BlockAvg = D2Traj;
      end
        
      % Sliding window processing
      if doSlidingWindowProcessing
        D2Traj = zeros(5,5,nSteps,nTraj);
        idx = 1:nSteps;
        for iTraj = 1:nTraj
          D2Traj(:,:,:,iTraj) = D2BlockAvg(:,:,idx,1);
          idx = idx + Par.lag;
          if alternatingTrajDirection
            idx = flip(idx);
          end
        end
      else
        D2Traj = D2BlockAvg;
      end
      
      % Cache processed trajectory if using a MD trajectory
      if useMD
        cacheD2Traj = D2Traj;
      end
    end
    
    % Combine local and global dynamics
    %---------------------------------------------------------------------------
    if ~isempty(Par.qLab)
      D2Lab = wigD(Par.qLab);
      D2Traj = multimatmult(D2Lab, D2Traj);
    end
        
    % Calculate Hamiltonians
    %---------------------------------------------------------------------------
    H0 = cacheTensors.H0;
    H = repmat(cacheTensors.Q0,[1,1,nSteps,nTraj]);
    for mp = 1:5
      for m = 1:5
        H = H + bsxfun(@times, D2Traj(m,mp,1:nSteps,:), cacheTensors.Q2{mp,m});
      end
    end
    
    % Calculate propagators
    %---------------------------------------------------------------------------
    algo = 'eig';
    useSzonly = true;
    if useSzonly
      % filter to eliminate Sx and Sz terms from Hamiltonians
      nStates = 2*Sys.Spins + 1;
      I(1:nStates(1)) = {ones(prod(nStates)/nStates(1))};
      SzFilter = blkdiag(I{:});
    else
      SzFilter = 1;
    end
    U0 = expm_(-1i*dtSpin*SzFilter.*H0,algo);  % interaction frame transformation operator
    U = zeros(size(H));
    for iStep = 1:nSteps
      for iTraj = 1:nTraj
        U_ = expm_(1i*dtSpin*SzFilter.*H(:,:,iStep,iTraj),algo);
        U(:,:,iStep,iTraj) = U0*U_;
      end
    end
    
    % Set up starting state of density matrix after pi/2 pulse, S_x
    %---------------------------------------------------------------------------
    rho = zeros(Sys.nStates,Sys.nStates,nSteps,nTraj);
    rho0 = sop(Sys.Spins,'x1');
    rho(:,:,1,:) = repmat(rho0,[1,1,1,nTraj]);
    
    % Propagate density matrix
    %---------------------------------------------------------------------------
    rho = propagate(rho,U,nSteps,PropagationMethod);
    
    % Average over trajectories (4th dimension)
    rho = mean(rho,4);
    
    % Multiply density matrix result by S_+ detection operator
    %---------------------------------------------------------------------------
    % Only keep the m_S=\beta subspace part that contributes to tr(S_{+}\rho(t))
    projector = sop(Sys.Spins,'+1');
    Sprho = multimatmult(repmat(projector,1,1,nSteps),rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise
    error('Quantum propagation method Opt.Method=''%s'' not recognized.',PropagationMethod);
end


end
%===============================================================================
%===============================================================================
%===============================================================================


%-------------------------------------------------------------------------------
% Propagate density matrix
function rho = propagate(rho,U,nSteps,Method)

switch Method
  
  case 'fast'
    % subspace propagation only requires U, but not U adjoint
    if size(rho,1)==1 % S=1/2 without nuclei
      for iStep = 2:nSteps
        rho(1,1,iStep,:) = U(1,1,iStep-1,:).^2.*rho(1,1,iStep-1,:);
      end
    else
      for iStep = 2:nSteps
        rho(:,:,iStep,:) = ...
          multimatmult( U(:,:,iStep-1,:),...
          multimatmult( rho(:,:,iStep-1,:),...
          U(:,:,iStep-1,:) ) );
      end
    end
    
  case 'ISTOs'
    
    Uadj = conj(permute(U,[2,1,3,4])); % adjoint of propagator
    for iStep = 2:nSteps
      rho(:,:,iStep,:) = ...
        multimatmult( U(:,:,iStep-1,:),...
        multimatmult( rho(:,:,iStep-1,:),...
        Uadj(:,:,iStep-1,:) ) );
    end
    
end

end

%-------------------------------------------------------------------------------
% calculate Wigner D-matrices of specified rank from quaternions for 
% rotation of ISTOs
function D2 = wigD(q)

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
