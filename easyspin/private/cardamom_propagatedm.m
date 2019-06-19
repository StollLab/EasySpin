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

persistent cacheTensors
persistent D2TrajMol
persistent gTensorState
persistent ATensorState
% persistent fullSteps
% persistent K2

% Preprocessing
% -------------------------------------------------------------------------

if ~isfield(Opt, 'Liouville')
  % Liouville space propagation is needed for the ensemble-averaged propagator
  Opt.Liouville = Opt.truncate;
end

Liouville = Opt.Liouville;
PropagationMethod = Opt.Method;

% Define a rotational dynamics time scale for integrating correlation functions
if ~isempty(MD)
  tauR = MD.tauR;
  isHMMfromMD = strcmp(Par.Model,'MD-HMM');
  useMD = true;
  isDirectfromMD = strcmp(Par.Model,'MD-direct');
else
  useMD = false;
  isHMMfromMD = false;
  if isfield(Sys,'Diff')       % TODO make this work for jumps and ISTOs
    tauR = 1/6/mean(Sys.Diff);
  end
end

if ~isfield(Par,'RTraj') && ~isfield(Par,'qTraj')
  error('Either Par.RTraj or Par.qTraj must be provided.')
end

Dt = Par.Dt;  % quantum propagation time step
nTraj = Par.nTraj;
doBlockAveraging = Par.isBlock;  % tensor time block averaging
nSteps = Par.nSteps;

truncate = Opt.truncate;  % ensemble averaged propagation

% store parameters for feeding into propagation function
Sim.nSteps = nSteps;
Sim.nTraj = nTraj;
Sim.truncate = truncate;
Sim.Liouville = Liouville;

% Simulation
% -------------------------------------------------------------------------

switch PropagationMethod
  case 'Nitroxide'  % see Ref [1]
    
    RTrajInv = permute(Par.RTraj,[2,1,3,4]);
    RLabInv = permute(Par.RLab,[2,1,3,4]);
    
    if ~isfield(Sys,'g')
      error('A g-tensor is required for the Nitroxide method.');
    end
    g = Sys.g;
    if ~isequal(size(g),[1,3])
      error('g-tensor must be a 3-vector.')
    end

    if ~isfield(Sys, 'A')
      error('An A-tensor is required for the Nitroxide method.');
    end
    A = Sys.A;
    if ~isequal(size(A),[1,3])
      error('A-tensor must be a 3-vector.')
    end
    
    if ~isHMMfromMD
      % Calculate time-dependent tensors from orientational trajectories
      gTensor = cardamom_tensortraj(g,Par.RTraj,RTrajInv);
      ATensor = cardamom_tensortraj(A,Par.RTraj,RTrajInv)*1e6*2*pi; % MHz (s^-1) -> Hz (rad s^-1)
      
    else
      % Calculate time-dependent tensors from state trajectories
      
      % Calculate the average interaction tensors for each state using
      % MD-derived frame trajectories and Viterbi state trajectories
      if isempty(gTensorState)
        % Perform MD-derived rotations on g- and A-tensors
        gTensorMD = cardamom_tensortraj(g,Par.RTraj,RTrajInv);
        ATensorMD = cardamom_tensortraj(A,Par.RTraj,RTrajInv)*1e6*2*pi; % MHz (s^-1) -> Hz (rad s^-1)
        
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
      
      for iState = 1:MD.nStates
        for iTraj = 1:nTraj
          idxState = Par.stateTraj(iTraj,:)==iState;
          gTensor(:,:,iTraj,idxState) = repmat(gTensorState(:,:,iState),[1,1,1,sum(idxState)]);
          ATensor(:,:,iTraj,idxState) = repmat(ATensorState(:,:,iState),[1,1,1,sum(idxState)]);
        end
      end
      
    end
    
    % Time block averaging and sliding window processing of tensors
    if doBlockAveraging
      gTensorBlock = zeros(3,3,size(gTensor,3),Par.nBlocks);
      ATensorBlock = zeros(3,3,size(ATensor,3),Par.nBlocks);

      % Average the interaction tensors over time blocks
      idx = 1:Par.BlockLength;
      for k = 1:Par.nBlocks
        gTensorBlock(:,:,:,k) = mean(gTensor(:,:,:,idx),4);
        ATensorBlock(:,:,:,k) = mean(ATensor(:,:,:,idx),4);
        idx = idx + Par.BlockLength;
      end

      if useMD && isDirectfromMD && nTraj>1
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
    if ~isempty(Par.RLab)
      gTensor = multimatmult(Par.RLab, multimatmult(gTensor, RLabInv));
      ATensor = multimatmult(Par.RLab, multimatmult(ATensor, RLabInv));
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
    Sprho = propagate(rho, U, PropagationMethod, Sim, Opt);

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
      
      if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end  % TODO include frames in cardamom
      
      [T,F,~,~,~] = magint(Sys,SpinOps,CenterField,false,false);
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
    
    % Process Wigner D-matrices
    % ---------------------------------------------------------------------
    % time block averaging and sliding window processing
    if doBlockAveraging
      if isempty(D2TrajMol)
        % if using an MD trajectory, it is best to process the molecular
        % dynamics once and store the result
        % this variable will be set to empty later if an MD trajectory is
        % not being used
        D2TrajMol = wigD(Par.qTraj);
        D2Avg = zeros(5,5,Par.nBlocks);
        
        % average the Wigner D-matrices over time blocks
        idx = 1:Par.BlockLength;
        for k = 1:Par.nBlocks
          D2Avg(:,:,k) = mean(D2TrajMol(:,:,:,idx),4);
          idx = idx + Par.BlockLength;
        end
        
        % perform sliding window processing if using MD trajectory explicitly
        if useMD
          if isDirectfromMD
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
      D2Traj = wigD(Par.qTraj);
    end
    
    % Check for new lab frame rotations
    if ~isempty(Par.qLab)
      D2Lab = wigD(Par.qLab);
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
        truncate = false;
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
      
      D2avg = mean(D2Traj, 4);  % time average of trajectories of Wigner matrices
      
    else
      % no approximation, use full propagator for all time steps
      fullSteps = nSteps;
      
    end
    
    % Calculate Hamiltonians
    % ---------------------------------------------------------------------
    H0 = cacheTensors.H0;
    H = repmat(cacheTensors.Q0,[1,1,nTraj,fullSteps]);
    for mp = 1:5
      for m = 1:5
        H = H + bsxfun(@times, D2Traj(m,mp,:,1:fullSteps), cacheTensors.Q2{mp,m});
      end
    end
    
    % Calculate propagators
    % ---------------------------------------------------------------------
    algo = 'eig';
    U0 = expm_(-1i*Dt*H0,algo);  % interaction frame transformation operator
    U = zeros(size(H));
    for iStep = 1:fullSteps
      for iTraj = 1:nTraj
        U(:,:,iTraj,iStep) = U0*expm_(1i*Dt*H(:,:,iTraj,iStep),algo);
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

    % Set up starting state of density matrix after pi/2 pulse, S_x
    % ---------------------------------------------------------------------
    Sim.fullSteps = fullSteps;
    rho = zeros(Sys.nStates,Sys.nStates,nTraj,nSteps);
    rho0 = sop(Sys.Spins,'x1');
    rho(:,:,:,1) = repmat(rho0,1,1,nTraj,1);
    
    % Propagate density matrix
    % ---------------------------------------------------------------------
    rho = propagate(rho, U, PropagationMethod, Sim, Opt);
    if Liouville
      if strcmp(Opt.debug.EqProp,'time')
        rho = reshape(rho,[Sys.nStates,Sys.nStates,nTraj,nSteps]);
      elseif strcmp(Opt.debug.EqProp,'all')
        rho = reshape(rho,[Sys.nStates,Sys.nStates,nSteps]);
      end
    end
    
    % Multiply density matrix result by S_+ detection operator
    % ---------------------------------------------------------------------
    % Only keep the m_S=\beta subspace part that contributes to tr(S_{+}\rho(t))
    projector = sop(Sys.Spins,'+1');    
    if strcmp(Opt.debug.EqProp,'time')
      Sprho = multimatmult(repmat(projector,1,1,nTraj,nSteps),rho);
    elseif strcmp(Opt.debug.EqProp,'all')
      Sprho = multimatmult(repmat(projector,1,1,nSteps),rho);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise
    error('Quantum propagation method Opt.Method=''%s'' not recognized.',PropagationMethod);
end


end
%===============================================================================
%===============================================================================
%===============================================================================

% Helper functions
% -------------------------------------------------------------------------

function rho = propagate(rho, U, Method, Sim, Opt)

fullSteps = Sim.fullSteps;
nSteps = Sim.nSteps;
nTraj = Sim.nTraj;
truncate = Sim.truncate;
%Liouville = Sim.Liouville;

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
    for iStep = 2:fullSteps
      rho(:,:,:,iStep) = ...
        multimatmult( U(:,:,:,iStep-1),...
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
      Uadj = conj(permute(U,[2,1,3,4])); % adjoint of propagator
      for iStep = 2:fullSteps
        rho(:,:,:,iStep) = ...
          multimatmult( U(:,:,:,iStep-1),...
          multimatmult( rho(:,:,:,iStep-1),...
          Uadj(:,:,:,iStep-1) ) );
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
