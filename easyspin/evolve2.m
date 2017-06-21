% evolve  Time domain evolution of density matrix
%
%   td = evolve(Sigma,Det,Ham,n,dt);
%   td = evolve(Sigma,Det,Ham,n,dt,IncScheme);
%   td = evolve(Sigma,Det,Ham,n,dt,IncScheme,Mix);
%
%   Evolves the density matrix Sigma under the Hamiltonian Ham with time
%   step dt n-1 times and detects using Det after each step. Hermitian
%   input matrices are assumed. td(1) is the value obtained by detecting
%   Sigma without evolution.
%
%   IncScheme determines the incrementation scheme and can be one of the
%   following (up to four sweep periods, up to two dimensions)
%
%     [1]           simple FID, 3p-ESEEM, echo transient, DEFENCE
%     [1 1]         2p-ESEEM, CP, 3p and 4p RIDME
%     [1 -1]        3p-DEER, 4p-DEER, PEANUT, 5p RIDME
%     [1 1 -1 -1]   SIFTER
%     [1 -1 -1 1]   7p-DEER
%     [1 -1 1 -1]
%
%     [1 2]         3p-ESEEM echo transient, HYSCORE, DONUT-HYSCORE
%     [1 2 1]       2D 3p-ESEEM
%     [1 1 2]       2p-ESEEM etc. with echo transient
%     [1 -1 2]      3p-DEER, 4p-DEER etc. with echo transient
%     [1 2 2 1]     2D CP
%     [1 2 -2 1]    2D PEANUT
%     [1 1 -1 -1 2] SIFTER with echo transient
%     [1 -1 -1 1 2] 7p-DEER with echo transient
%
%   [1] is the default. For an explanation of the format, see the
%   documentation.
%
%   Mix is a cell array containing the propagators of the mixing
%   block(s), for experiments with more than 1 sweep period.
%
%   td is a vector/matrix of the signal with t1 along dimension 1
%   and t2 along dimension 2.



function [TimeArray, SignalArray, Sigma, DensityMatrices, Events] = evolve2(Sigma,Ham0,Det,Events,Relaxation,Vary)

if (nargin==0), help(mfilename); return; end

% if (nargout<1), error('Not enough output arguments.'); end
if (nargout>5), error('Too many output arguments.'); end
if (nargin<4) || (nargin>6), error('Wrong number of input arguments!'); end

%move this this into Events?
% if (nargin<6), IncScheme = 1; end
% if (nargin<7), Mix = {}; end
%
% if any(mod(n,1)~=0) || any(n<=0)
%     error('n, the number of points (4th argument), must be a positive integer.');
% end

method = 'stepwise';

nEvents = length(Events);

% questions:
% - we have to look into the normalization for the Detection operators. +
% or - Detection operators need to divided by two before normalizations are
% computed. what should the normalization be for transition selective
% operators or for detecting populations
% - If no detection is requested, I can combine all propagators to one single
% one: Utotal = U1*U2*U3.... can this also be done for Liouvillians

switch method
  
  case 'stepwise'
    
    nDet = numel(Det);
    normsDet = zeros(1,nDet);
    Ham0 = Ham0*2*pi;
    
    [nAcquisition, nDimensions] = size(Vary.Table);
    
    initialSigma = Sigma;
    
    %%%% This checks if all time traces will have the same length and can 
    % therefore be stored in an array
    StoreInArray = false;
    for iDimension = 1 : nDimensions
      for Event2Check = Vary.Events{iDimension}
        if Events{Event2Check}.Detection == 1
          nDataSets = length(Vary.Dimension{iDimension}.ts{Event2Check});
          sizetocompare = size(Vary.Dimension{iDimension}.ts{Event2Check}{1},2);
          for i = 2 : nDataSets
            if sizetocompare ~= size(Vary.Dimension{iDimension}.ts{Event2Check}{i},2)
              StoreInArray = true;
            end
          end
        end
      end
    end
    
    if StoreInArray
      SignalArray = cell(1,nAcquisition);
      TimeArray = SignalArray;
    end
    %----------------------------------------------------------------------
    
    for iAcquisition = 1 : nAcquisition
      tic
      
      Sigma = initialSigma;
      
      for iDimension = 1 : nDimensions
        iPoint = Vary.Table(iAcquisition,iDimension);
        
        for iEvent = Vary.Events{iDimension}
                   
          switch Events{iEvent}.type
            case 'pulse'
              Events{iEvent}.IQ = Vary.Dimension{iDimension}.IQs{iEvent}{iPoint};
              Events{iEvent}.t = Vary.Dimension{iDimension}.ts{iEvent}{iPoint};
            case 'free evolution'
              Events{iEvent}.t = Vary.Dimension{iDimension}.ts{iEvent}{iPoint};
          end
          
          Events{iEvent}.Propagation = [];
          
          
        end
        
        
      end
      ttotal = 0;
      firstDetection = 1;
      startTrace = 2;
      
      
      
      
      for iEvent = 1 : nEvents
        currentEvent = Events{iEvent};
        
        if isempty(currentEvent.Propagation) || (~isfield(currentEvent.Propagation,'Utotal') || ~isfield(currentEvent.Propagation,'Ltotal'))
          
          
          switch currentEvent.type
            case 'pulse'
              
              dt = currentEvent.t(2) - currentEvent.t(1);
              
              nPhaseCycle = size(currentEvent.PhaseCycle,1);
              %----------------------------------------------------------------
              % convert the IQ wave form into a binary form, so that it
              % is possible to use a propagator look up table
              %----------------------------------------------------------------
              
              rMax = max(abs(real(currentEvent.IQ(:)))); % max(abs(real....(:)))
              iMax = max(abs(imag(currentEvent.IQ(:))));
              
              if rMax > iMax
                MaxWave = rMax;
              else
                MaxWave = iMax;
              end
              
              vertRes = 1024;
              
              scale = 2*2*2*pi*MaxWave/vertRes;
              % one factor 2 is required because of linearly polarized irradiation
              % the other because of the digitization of the wave
              
              tvector(1) = 0;
              tvector(2:length(currentEvent.t)+1) = currentEvent.t+dt;
              
              realArbitrary = real(currentEvent.IQ)/MaxWave;
              realBinary = floor(vertRes*(realArbitrary+1)/2);
              realBinary(realBinary == vertRes) = vertRes-1;
              
              if currentEvent.ComplexExcitation == 1
                imagArbitrary = imag(currentEvent.IQ)/MaxWave;
                imagBinary = floor(vertRes*(imagArbitrary+1)/2);
                imagBinary(imagBinary==vertRes) = vertRes-1;
              end
              
              if ~currentEvent.Relaxation
                
                % Calculation or Loading, if possible, of Propagators
                if currentEvent.ComplexExcitation == 0
                  if ~isempty(currentEvent.Propagation) && isfield(currentEvent.Propagation,'UTable')
                    UTable = currentEvent.Propagation.UTable;
                  else
                    
                    UTable = buildPropagators(Ham0,currentEvent.xOp,dt,vertRes,scale);
                    
                    Events{iEvent}.Propagation.UTable = UTable;
                  end
                  
                end
                
                if currentEvent.Detection == 1
                  Utotal = cell(nPhaseCycle,length(realBinary));
                end
                
                % Loops over the Phasecycles
                
                for iPhaseCycle = 1 : nPhaseCycle
                  % Propagation for one waveform
                  for iWavePoint = 1 : length(realBinary)
                    if currentEvent.ComplexExcitation == 0
                      % Load propagators if Complex Excitation is off
                      U = UTable{realBinary(iPhaseCycle,iWavePoint)+1};
                    else % For active Complex Excitation Propagators need to be recalculated
                      Ham1 = scale*(realBinary(iPhaseCycle,iWavePoint)-vertRes/2)*real(currentEvent.xOp)+scale*(imagBinary(iPhaseCycle,iWavePoint)-vertRes/2)*imag(currentEvent.xOp);
                      Ham =  Ham0+Ham1;
                      U = Propagator(Ham,dt);
                    end
                    
                    if currentEvent.Detection == 1
                      Utotal{iPhaseCycle,iWavePoint} = U;
                    else
                      if iWavePoint == 1
                        Utotal{iPhaseCycle} = U;
                      else
                        Utotal{iPhaseCycle} = U*Utotal{iPhaseCycle};
                      end
                    end
                  end
                end
                
                Events{iEvent}.Propagation.Utotal = Utotal;
                
              elseif currentEvent.Relaxation
                
                n = size(Sigma,1);
                Gamma = Relaxation.Gamma;
                equilibriumState = reshape(Relaxation.equilibriumState,n*n,1);
                
                if ~isfield(currentEvent.Propagation,'Ltotal')
                  % Calculation or Loading, if possible, of Propagators
                  if currentEvent.ComplexExcitation == 0
                    if ~isempty(currentEvent.Propagation) && isfield(currentEvent.Propagation,'LTable')
                      LTable = currentEvent.Propagation.LTable;
                      SigmassTable = currentEvent.Propagation.SigmassTable;
                    else
                      
                      [LTable, SigmassTable] = buildLiouvillians(Ham0,Gamma,equilibriumState,currentEvent.xOp,dt,vertRes,scale);
                      
                      Events{iEvent}.Propagation.LTable = LTable;
                      Events{iEvent}.Propagation.SigmassTable = SigmassTable;
                    end
                    
                  end
                  
                  if currentEvent.Detection == 1
                    Ltotal = cell(nPhaseCycle,length(realBinary));
                    SigmaSStotal = Ltotal;
                  end
                  
                  % Loops over the Phasecycles
                  
                  for iPhaseCycle = 1 : nPhaseCycle
                    % Propagation for one waveform
                    for iWavePoint = 1 : length(realBinary)
                      if currentEvent.ComplexExcitation == 0
                        % Load Liouvillians if Complex Excitation is off
                        L=LTable{realBinary(iPhaseCycle,iWavePoint)+1};
                        SigmaSS=SigmassTable{realBinary(iPhaseCycle,iWavePoint)+1};
                      else
                        % if complex excitation is requested, usage of tables
                        % is not feasible, and Liouvillians and state state
                        % density matrices are computed for each time step Ham1 = scale*(realBinary(iPhaseCycle,iWavePoint)-vertRes/2)*real(currentEvent.xOp)+scale*(imagBinary(iPhaseCycle,iWavePoint)-vertRes/2)*imag(currentEvent.xOp);
                        Ham1 = scale*(realBinary(iPhaseCycle,iWavePoint)-vertRes/2)*real(currentEvent.xOp)+scale*(imagBinary(iPhaseCycle,iWavePoint)-vertRes/2)*imag(currentEvent.xOp);
                        Ham = Ham0+Ham1;
                        [L, SigmaSS] = Liouvillian(Ham,Gamma,equilibriumState,dt);
                      end
                      
                      if currentEvent.Detection == 1
                        Ltotal{iPhaseCycle,iWavePoint} = L;
                        SigmaSStotal{iPhaseCycle,iWavePoint} = SigmaSS;
                      else
                        Ltotal{iPhaseCycle,iWavePoint} = L;
                        SigmaSStotal{iPhaseCycle,iWavePoint} = SigmaSS;
                        %                   Here the total Liouvillian and SigmaSS need to be
                        %                   calculated, if possibly
                      end
                    end
                  end
                  
                  Events{iEvent}.Propagation.Ltotal = Ltotal;
                  Events{iEvent}.Propagation.SigmaSStotal = SigmaSStotal;
                end
              end
              
          end
          
        end
      end
      
      for iEvent = 1 : nEvents
        
        currentEvent = Events{iEvent};
        
        if length(currentEvent.t) == 1
          dt = currentEvent.t;
          tvector = [0 currentEvent.t];
        else
          dt = currentEvent.t(2) - currentEvent.t(1);
          tvector = currentEvent.t;
        end
        
        if currentEvent.Detection == 1
          currentSignal = zeros(nDet,length(tvector));
          n = size(Sigma,1);
          
          for iDet = 1:length(Det)
            
            Det{iDet} = reshape(Det{iDet}.',1,n^2);
            normsDet(iDet) = Det{iDet}*Det{iDet}';
            %           normsDet(kk) = 1;
            Density = Sigma(:);
            %           normsDet(kk) = sum(sum(Det{kk}.*Det{kk}));
            %           currentSignal(kk,1) = sum(sum(Det{kk}.*Sigma.'))/normsDet(kk);
            currentSignal(iDet,1) = Det{iDet}*Density/normsDet(iDet);
          end
          
        else
          currentSignal=[];
        end
        
        if currentEvent.storeDensityMatrix == 1
          DensityMatrices=cell(1,length(currentEvent.t));
          DensityMatrices{1}=Sigma;
        else
          DensityMatrices = [];
        end
        
        
        switch currentEvent.type
          case 'pulse'
            
            tvector(1) = 0;
            tvector(2:length(currentEvent.t)+1) = currentEvent.t+dt;
            
            nPhaseCycle = size(currentEvent.PhaseCycle,1);
            
            if nPhaseCycle>1
              StateBeforePC = Sigma;
              loopState = zeros(size(Sigma));
              PCnorm = sum(abs(currentEvent.PhaseCycle(:,2)));
            end
            
            %----------------------------------------------------------------
            % Propagation Starts Here
            %----------------------------------------------------------------
            
            if ~currentEvent.Relaxation
              
              
              % Loops over the Phasecycles
              for iPhaseCycle = 1 : nPhaseCycle
                
                
                if currentEvent.Detection == 0
                  U = currentEvent.Propagation.Utotal{iPhaseCycle};
                  Sigma = U*Sigma*U';
                else
                  for iWavePoint = 1 : size(currentEvent.Propagation.Utotal,2)
                    U = currentEvent.Propagation.Utotal{iPhaseCycle,iWavePoint};
                    Sigma = U*Sigma*U';
                    
                    % Computes Expectation Values if requested
                    if currentEvent.Detection == 1
                      for iDet = 1:nDet
                        Density = Sigma(:);
                        currentSignal(iDet,iWavePoint+1) = Det{iDet}*Density/normsDet(iDet);
                        %                   currentSignal(j,k+1)=sum(sum(Det{j}.*Sigma.'))/normsDet(j);
                      end
                    end
                    
                    % Store Density Matrices if requested
                    if currentEvent.storeDensityMatrix == 1
                      DensityMatrices{iWavePoint+1}=Sigma;
                    end
                    
                  end
                end
                
                % Combine Results from current phase cycle with previous ones
                if nPhaseCycle > 1
                  PCweight = currentEvent.PhaseCycle(iPhaseCycle,2)/PCnorm;
                  loopState = loopState+PCweight*Sigma;
                  if currentEvent.Detection == 1
                    if iPhaseCycle ~= 1
                      PCSignal = PCSignal+PCweight*currentSignal;
                    else
                      PCSignal = PCweight*currentSignal;
                    end
                  end
                  if iPhaseCycle ~= nPhaseCycle
                    Sigma = StateBeforePC;
                  end
                end
                
              end
              
            elseif currentEvent.Relaxation % Propagation in Liouville space
              
              
              % Loops over the Phasecycles
              for iPhaseCycle = 1 : nPhaseCycle
                SigmaVector = reshape(Sigma,n*n,1);
                
                if currentEvent.Detection == 3 % if propagation in Liouville Space in one step is possible, this has to go here
                  Utotal = currentEvent.Propagation.Utotal{iPhaseCycle};
                  Sigma = Utotal*Sigma*Utotal';
                else
                  for iWavePoint = 1 : size(currentEvent.Propagation.Ltotal,2)
                    L = currentEvent.Propagation.Ltotal{iPhaseCycle,iWavePoint};
                    SigmaSS = currentEvent.Propagation.SigmaSStotal{iPhaseCycle,iWavePoint};
                    
                    SigmaVector = SigmaSS+L*(SigmaVector-SigmaSS);
                    
                    
                    % Computes Expectation Values if requested
                    if currentEvent.Detection == 1
                      Sigma = reshape(SigmaVector,n,n);
                      for iDet = 1:nDet
                        Density = Sigma(:);
                        currentSignal(iDet,iWavePoint+1) = Det{iDet}*Density/normsDet(iDet);
                        %                   currentSignal(j,k+1)=sum(sum(Det{j}.*Sigma.'))/normsDet(j);
                      end
                    end
                    
                    % Store Density Matrices if requested
                    if currentEvent.storeDensityMatrix == 1
                      Sigma = reshape(SigmaVector,n,n);
                      DensityMatrices{iWavePoint+1}=Sigma;
                    end
                    
                  end
                end
                
                % Combine Results from current phase cycle with previous ones
                if nPhaseCycle > 1
                  PCweight = currentEvent.PhaseCycle(iPhaseCycle,2)/PCnorm;
                  loopState = loopState+PCweight*Sigma;
                  if currentEvent.Detection == 1
                    if iPhaseCycle ~= 1
                      PCSignal = PCSignal+PCweight*currentSignal;
                    else
                      PCSignal = PCweight*currentSignal;
                    end
                  end
                  if iPhaseCycle ~= nPhaseCycle
                    Sigma = StateBeforePC;
                  end
                end
                
              end
              
            end
            
            % After Propagation and if PhaseCycling was active, the results
            % from phase cycling are returned
            if nPhaseCycle > 1
              Sigma = loopState;
              if currentEvent.Detection == 1
                currentSignal = PCSignal;
              end
            end
            
          case 'free evolution'
            
            if ~currentEvent.Relaxation
              % If Detection is off during evolution, the entire evolution
              % can be propagated in one step
              if currentEvent.Detection == 0
                dt = currentEvent.t(end);
                tvector = [0 dt];
              end
              
              if isfield(currentEvent.Propagation,'Utotal')
                U = currentEvent.Propagation.Utotal;
              else
                U = Propagator(Ham0,dt);
                Events{iEvent}.Propagation.Utotal = U;
              end
              
              % Propagation starts here
              for itvector=2:length(tvector)
                
                Sigma=U*Sigma*U';
                
                if currentEvent.Detection == 1
                  for iDet = 1:nDet
                    Density = Sigma(:);
                    currentSignal(iDet,itvector) = Det{iDet}*Density/normsDet(iDet);
                    %               currentSignal(j,k)=sum(sum(Det{j}.*Sigma.'))/normsDet(j);
                  end
                end
                
                if currentEvent.storeDensityMatrix == 1
                  DensityMatrices{itvector}=Sigma;
                end
                
              end
              
            elseif currentEvent.Relaxation % Propagate in Liouville space
              
              n = size(Sigma,1);
              Gamma = Relaxation.Gamma;
              equilibriumState = reshape(Relaxation.equilibriumState,n*n,1);
              SigmaVector = reshape(Sigma,n*n,1);
              
              if currentEvent.Detection == 0
                dt = currentEvent.t(end);
                tvector = [0 dt];
              end
              
              if isfield(currentEvent.Propagation,'Ltotal')
                L = currentEvent.Propagation.Ltotal;
                SigmaSS = currentEvent.Propagation.SigmaSStotal;
              else
                [L, SigmaSS] = Liouvillian(Ham0,Gamma,equilibriumState,dt);
                Events{iEvent}.Propagation.Ltotal = L;
                Events{iEvent}.Propagation.SigmaSStotal = SigmaSS;
              end
              
              
              % Propagation
              for itvector = 2:length(tvector)
                SigmaVector = SigmaSS+L*(SigmaVector-SigmaSS);
                Sigma = reshape(SigmaVector,n,n);
                
                if currentEvent.Detection == 1
                  for iDet = 1:nDet
                    Density = Sigma(:);
                    currentSignal(iDet,itvector) = Det{iDet}*Density/normsDet(iDet);
                    %                 currentSignal(j,k)=sum(sum(Det{j}.*Sigma.'))/normsDet(j);
                  end
                end
                
                if currentEvent.storeDensityMatrix == 1
                  DensityMatrices{itvector} = Sigma;
                end
              end
              
              
            end
            
        end
        
        % This combines the signals and time axes from all detected event
        if firstDetection && ~isempty(currentSignal)
          Signal = [];
          t = [];
          % store first point of the first signal
          Signal(:,1) = currentSignal(:,1);
          t(1,1) = ttotal;
          
          % now add all the others timepoints
          nSignal = size(currentSignal,2);
          endTrace = startTrace+nSignal-2;
          Signal(:,startTrace:endTrace) = currentSignal(:,2:end);
          t(1,startTrace:endTrace) = tvector(2:end)+ttotal;
          startTrace = endTrace+1;
          firstDetection = 0;
          
        elseif ~isempty(currentSignal)
          % adding other signals, this confusing index is necessary in
          % order to avoid double counting last point of a signal and the
          % first point of the succiding signal
          nSignal = size(currentSignal,2);
          endTrace = startTrace+nSignal-2;
          Signal(:,startTrace:endTrace) = currentSignal(:,2:end);
          t(1,startTrace:endTrace) = tvector(2:end)+ttotal;
          startTrace = endTrace+1;
        end
        
        % Update Total Time, necessary to keep correct timings  if events are
        % not detected
        ttotal = ttotal + tvector(end);
        
      end
      
      if ~firstDetection
        % move all the detected signals into large array here!
        if StoreInArray
          SignalArray{iAcquisition} = Signal;
          TimeArray{iAcquisition} = t;
        else
          if iAcquisition == 1
            SignalSize = size(Signal);
            SignalArray = zeros(SignalSize(1),SignalSize(2),nAcquisition);
            SignalArray(:,:,iAcquisition) = Signal;
            TimeArray = t;
          else
            SignalArray(:,:,iAcquisition) = Signal;
          end
        end
      end
      toc
    end
    
    if firstDetection
      SignalArray = [];
      TimeArray = [];
    end
    

  case 'incrementation scheme'
    
    
    
    % IncScheme check
    %------------------------------------------------------------
    if (length(IncScheme)>1) && (nargin<7),
      error('The requested IncScheme requires mixing propagators, but none are provided!');
    end
    if any((abs(IncScheme)~=1) & (abs(IncScheme)~=2))
      error('IncScheme can contain only 1, -1, 2, and -2.');
    end
    
    nEvolutionPeriods = numel(IncScheme);
    nDimensions = max(abs(IncScheme));
    
    % Parameter parsing
    %------------------------------------------------------------
    if ~iscell(Det)
      Det = {Det};
    end
    nDetectors = numel(Det);
    
    if ~iscell(Mix)
      Mix = {Mix};
    end
    nMixingBlocks = numel(Mix);
    
    if (nMixingBlocks~=nEvolutionPeriods-1),
      error('Number of mixing propagators not correct! %d are needed.',nEvolutionPeriods-1);
    end
    N = size(Sigma,1);
    
    if (nDimensions==1)
      for iDet = 1:nDetectors
        Signal{iDet} = zeros(n,1);
      end
      if iscell(Ham0), Ham0 = Ham0{1}; end
    else
      if numel(dt)==1, dt = [dt dt]; end
      if numel(n)==1, n = [n n]; end
      for iDet = 1:nDetectors
        Signal{iDet} = zeros(n);
      end
    end
    if nDetectors==1
      Signal = Signal{1};
    end
    
    % Transform all operators to Hamiltonian eigenbasis (eigenbases)
    %---------------------------------------------------------------
    if ~iscell(Ham0)
      
      if nnz(Ham0)==nnz(diag(Ham0)) % Check if Hamiltonian is already diagonal
        E = diag(Ham0);
        Density = Sigma;
        Detector = Det;
      else
        % Diagonalize Hamiltonian
        [Vecs,E] = eig(Ham0); % MHz, E doesn't have to be sorted
        E = real(diag(E));
        % Transform all other matrices to Hamiltonian eigenbasis
        Density = Vecs'*Sigma*Vecs;
        for iMix = 1:nMixingBlocks
          Mix{iMix} = Vecs'*Mix{iMix}*Vecs;
        end
        for iDet = 1:nDetectors
          Detector{iDet} = Vecs'*Det{iDet}*Vecs;
        end
      end
      
      % Define free evolution propagators
      if (nDimensions==1)
        diagU = exp(-2i*pi*dt*E);
      else
        diagUX = exp(-2i*pi*dt(1)*E);
        diagUY = exp(-2i*pi*dt(2)*E);
      end
      
    else
      
      % Check if Hamiltonians are already diagonal
      if (nnz(Ham0{1})==nnz(diag(Ham0{1})) && nnz(Ham0{1})==nnz(diag(Ham0{1})))
        Ex = Ham0{1};
        Ey = Ham0{2};
        Density = Sigma;
        Detector = Det;
      else
        
        % Diagonalize Hamiltonians
        [Vecs{1},Ex] = eig(Ham0{1});
        [Vecs{2},Ey] = eig(Ham0{2});
        
        % Transform all other matrices to Hamiltonian eigenbasis
        d = abs(IncScheme);
        Density = Vecs{d(1)}'*Sigma*Vecs{d(1)};
        for iMix = 1:nMixingBlocks
          Mix{iMix} = Vecs{d(iMix+1)}'*Mix{iMix}*Vecs{d(iMix)};
        end
        for iDet = 1:nDetectors
          Detector{iDet} = Vecs{d(end)}'*Det{iDet}*Vecs{d(end)};
        end
        
      end
      
      % Define free evolution propagators
      diagUX = exp((-2i*pi*dt(1))*real(diag(Ex)));
      diagUY = exp((-2i*pi*dt(2))*real(diag(Ey)));
      
    end
    
    % Time-domain evolution, IncScheme switchyard
    %------------------------------------------------------------
    % The following implementations in the propagator eigenframes
    % (giving diagonal propagators) make use of the following simplifications
    % of the matrix multiplications associated with the propagations:
    %
    %   U*Density*U'   = (diagU*diagU').*Density
    %   U*Propagator*U = (diagU*diagU.').*Propagator
    %   (U^-1)*Propagator*U = U'*Propagator*U = (conj(diagU)*diagU.').*Propagator
    %   U*Propagator*(U^-1) = U*Propagator*U' = (diagU*diagU').*Propagator
    %
    % where U are diagonal matrices and diagU are vectors of eigenvalues.
    
    % Pre-reshape for trace calculation
    for iDet = 1:nDetectors
      Detector{iDet} = reshape(Detector{iDet}.',1,N^2);
    end
    if nDetectors==1
      Detector = Detector{1};
    end
    
    if isequal(IncScheme,1) % IncScheme [1]
      FinalDensity = Density(:);
      U_ = diagU*diagU';
      U_ = U_(:);
      for ix = 1:n
        % Compute trace(Detector*FinalDensity)
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity;
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity;
          end
        end
        FinalDensity = U_.*FinalDensity; % equivalent to U*FinalDensity*U'
      end
      
    elseif isequal(IncScheme,[1 1]) % IncScheme [1 1]
      UU_ = diagU*diagU.';
      % It is not necessary to evolve the initial density matrix.
      % Only the mixing propagator needs to be evolved.
      Mix1 = Mix{1};
      for ix = 1:n
        % Compute density right before Detection
        FinalDensity = Mix1*Density*Mix1';
        % Compute trace(Detector*FinalDensity)
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity(:);
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
          end
        end
        Mix1 = UU_.*Mix1; % equivalent to U*Mix1*U
      end
      
    elseif isequal(IncScheme,[1 -1]) % IncScheme [1 -1]
      %   % Pre-propagate mixing propagator to end of second period (= start of
      %   % experiment)
      %   MixX = diag(diagU.^n)*Mix{1};
      MixX =Mix{1};
      UtU_ = conj(diagU)*diagU.';
      for ix = 1:n
        FinalDensity = MixX*Density*MixX';
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity(:);
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
          end
        end
        MixX = UtU_.*MixX; % equivalent to U^-1*MixX*U
      end
      
    elseif isequal(IncScheme,[1 2]) % IncScheme [1 2]
      UX_ = diagUX*diagUX';
      UY_ = diagUY*diagUY';
      UY_ = reshape(UY_,N^2,1);
      Mix1 = Mix{1};
      for ix = 1:n(1)
        FinalDensity = reshape(Mix1*Density*Mix1',N^2,1);
        for iy = 1:n(2)
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity;
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity;
            end
          end
          FinalDensity = UY_.*FinalDensity; % equivalent to UY*Density*UY';
        end
        Density = UX_.*Density; % equivalent to UX*Density*UX';
      end
      
    elseif isequal(IncScheme,[1 1 2]) % IncScheme [1 1 2]
      UUX_ = diagUX*diagUX.';
      UY_ = diagUY*diagUY';
      Mix1 = Mix{1};
      Mix2 = Mix{2};
      for ix = 1:n(1)
        M = Mix2*Mix1;
        FinalDensity = M*Density*M';
        for iy = 1:n(2)
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
        end
        Mix1 = UUX_.*Mix1; % equivalent to UX*Mix1*UX
      end
      
    elseif isequal(IncScheme,[1 -1 2]) % IncScheme [1 -1 2]
      UtUX_ = conj(diagUX)*diagUX.';
      UY_ = diagUY*diagUY';
      %   % Pre-propagate mixing propagator to end of second period (= start of
      %   % experiment)
      %   MixX = diag(diagUX.^n(1))*Mix{1};
      MixX = Mix{1};
      Mix2 = Mix{2};
      for ix = 1:n(1)
        M = Mix2*MixX;
        FinalDensity = M*Density*M';
        for iy = 1:n(2)
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
        end
        MixX = UtUX_.*MixX; % equivalent to U^-1*MixX*U
      end
      
    elseif isequal(IncScheme,[1 2 1]) % IncScheme [1 2 1]
      Mix1 = Mix{1};
      Mix2 = Mix{2};
      UUX_ = diagUX*diagUX.';
      UY = diag(diagUY);
      for iy = 1:n(2)
        MixY = Mix2*Mix1;
        MixYadj = MixY';
        for ix = 1:n(1)
          FinalDensity = MixY*Density*MixYadj;
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          MixY = UUX_.*MixY; % equivalent to UX*MixY*UX
        end
        Mix1 = UY*Mix1;
      end
      
    elseif isequal(IncScheme,[1 2 2 1]) % IncScheme [1 2 2 1]
      Mix1 = Mix{1};
      Mix2 = Mix{2};
      Mix3 = Mix{3};
      UUX_ = diagUX*diagUX.';
      UUY_ = diagUY*diagUY.';
      for iy = 1:n(2)
        MixY = Mix3*Mix2*Mix1;
        MixYadj = MixY';
        for ix = 1:n(1)
          FinalDensity = MixY*Density*MixYadj;
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          MixY = UUX_.*MixY; % equivalent to UX*MixY*UX
        end
        Mix2 = UUY_.*Mix2; % equivalent to UY*Mix2*UY
      end
      
    elseif isequal(IncScheme,[1 2 -2 1]) % IncScheme [1 2 -2 1]
      Mix1 = Mix{1};
      %   Mix2 = diag(diagUY.^n(2))*Mix{2}; % pre-propagate to endpoint of third delay
      Mix2 = Mix{2};
      Mix3 = Mix{3};
      UUX_ = diagUX*diagUX.';
      UtUY_ = conj(diagUY)*diagUY.';
      for iy = 1:n(2)
        MixY = Mix3*Mix2*Mix1;
        MixYadj = MixY';
        for ix = 1:n(1)
          FinalDensity = MixY*Density*MixYadj;
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          MixY = UUX_.*MixY; % equivalent to UX*MixY*UX
        end
        Mix2 = UtUY_.*Mix2; % equivalent to UY'*Mix2*UY
      end
      
    elseif isequal(IncScheme,[1 -1 1 -1]) % IncScheme [1 -1 1 -1]
      %   Mix1X = diag(diagU.^n)*Mix{1}; % pre-propagate to endpoint of second delay
      Mix1X = Mix{1};
      Mix2 = Mix{2};
      %   Mix3X = diag(diagU.^n)*Mix{3}; % pre-propagate to endpoint of fourth delay
      Mix3X = Mix{3};
      UtU_ = conj(diagU)*diagU.'; % propagator for Mix1 and Mix3 (add before, remove after)
      for ix = 1:n
        MixX = Mix3X*Mix2*Mix1X;
        FinalDensity = MixX*Density*MixX';
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity(:);
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
          end
        end
        Mix1X = UtU_.*Mix1X; % equivalent to U'*Mix1X*U
        Mix3X = UtU_.*Mix3X; % equivalent to U'*Mix3X*U
      end
      
    elseif isequal(IncScheme,[1 1 -1 -1]) % IncScheme [1 1 -1 -1]
      Mix1 = Mix{1};
      %   Mix2X = diag(diagU.^n)*Mix{2}; % pre-propagate to endpoint of third delay
      %   Mix3X = diag(diagU.^n)*Mix{3}; % pre-propagate to endpoint of fourth delay
      Mix2X = Mix{2};
      Mix3X = Mix{3};
      UU1_ = diagU*diagU.'; % propagator for Mix1 (add before and after)
      UU3_ = conj(diagU*diagU.'); % propagator for Mix3 (remove before and after)
      for ix = 1:n
        MixX = Mix3X*Mix2X*Mix1;
        FinalDensity = MixX*Density*MixX';
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity(:);
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
          end
        end
        Mix1 = UU1_.*Mix1; % equivalent to U*Mix1*U
        Mix3X = UU3_.*Mix3X; % equivalent to U'*Mix3X*U'
      end
      
    elseif isequal(IncScheme,[1 1 -1 -1 2]) % IncScheme [1 1 -1 -1 2]
      Mix1 = Mix{1};
      %   Mix2X = diag(diagUX.^n(1))*Mix{2}; % pre-propagate to endpoint of third delay
      %   Mix3X = diag(diagUX.^n(1))*Mix{3}; % pre-propagate to endpoint of fourth delay
      Mix2X = Mix{2};
      Mix3X = Mix{3};
      Mix4 = Mix{4};
      UU1_ = diagUX*diagUX.'; % propagator for Mix1 (add before and after)
      UU3_ = conj(diagUX*diagUX.'); % propagator for Mix3 (remove before and after)
      UY_ = diagUY*diagUY';
      for ix = 1:n(1)
        MixX = Mix4*Mix3X*Mix2X*Mix1;
        FinalDensity = MixX*Density*MixX';
        for iy = 1:n(2)
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
        end
        Mix1 = UU1_.*Mix1; % equivalent to U*Mix1*U
        Mix3X = UU3_.*Mix3X; % equivalent to U'*Mix3X*U'
      end
      
    elseif isequal(IncScheme,[1 -1 -1 1]) % IncScheme [1 -1 -1 1]
      %   Mix1X = diag(diagU.^n)*Mix{1}; % pre-propagate to endpoint of second delay
      Mix1X = Mix{1};
      Mix2 = Mix{2};
      %   Mix3X = Mix{3}*diag(diagU.^n); % forward-propagate to endpoint of fourth delay
      Mix3X = Mix{3};
      UtU1_ = conj(diagU)*diagU.'; % propagator for Mix1 (add before, remove after)
      UtU3_ = diagU*diagU'; % propagator for Mix3 (remove before, add after)
      for ix = 1:n
        MixX = Mix3X*Mix2*Mix1X;
        FinalDensity = MixX*Density*MixX';
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity(:);
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
          end
        end
        Mix1X = UtU1_.*Mix1X; % equivalent to U'*Mix1X*U
        Mix3X = UtU3_.*Mix3X; % equivalent to U*Mix3X*U'
      end
      
    elseif isequal(IncScheme,[1 -1 -1 1 2]) % IncScheme [1 -1 -1 1 2]
      %   Mix1X = diag(diagUX.^n(1))*Mix{1}; % pre-propagate to endpoint of second delay
      Mix1X = Mix{1};
      Mix2 = Mix{2};
      %   Mix3X = Mix{3}*diag(diagUX.^n(1)); % forward-propagate to endpoint of fourth delay
      Mix3X = Mix{3};
      Mix4 = Mix{4};
      UtU1_ = conj(diagUX)*diagUX.'; % propagator for Mix1 (add before, remove after)
      UtU3_ = diagUX*diagUX'; % propagator for Mix3 (remove before, add after)
      UY_ = diagUY*diagUY';
      for ix = 1:n(1)
        MixX = Mix4*Mix3X*Mix2*Mix1X;
        FinalDensity = MixX*Density*MixX';
        for iy = 1:n(2)
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
        end
        Mix1X = UtU1_.*Mix1X; % equivalent to U'*Mix1X*U
        Mix3X = UtU3_.*Mix3X; % equivalent to U*Mix3X*U'
      end
      
    else
      error('Unsupported incrementation scheme!');
    end
    
    return
end


function UTable = buildPropagators(Ham0,xOp,dt,vertRes,scale)

UTable = cell(1,vertRes);

for iRes = 0:vertRes-1
  Ham1 = scale*(iRes-vertRes/2)*xOp;
  Ham = Ham0+Ham1;
  U = Propagator(Ham,dt);
  UTable{iRes+1} = U;
end

function U = Propagator(Ham,dt)
U = expm_fastc(-1i*Ham*dt);

function [LTable, SigmassTable] = buildLiouvillians(Ham0,Gamma,equilibriumState,xOp,dt,vertRes,scale)

SigmassTable = cell(1,vertRes);
LTable = cell(1,vertRes);

for ivertRes = 0:vertRes-1
  Ham1 = scale*(ivertRes-vertRes/2)*xOp;
  Ham = Ham0 + Ham1;
  [L, SigmaSS] = Liouvillian(Ham,Gamma,equilibriumState,dt);
  LTable{ivertRes+1} = L;
  SigmassTable{ivertRes+1} = SigmaSS;
end

function [L, SigmaSS] = Liouvillian(Ham,Gamma,equilibriumState,dt)

n = size(Ham,1);

HamSuOp = kron(eye(n,n),Ham)-kron(Ham.',eye(n,n));
L = -1i*HamSuOp-Gamma;
SigmaSS = Gamma*equilibriumState; % steady state solutions for the density matrices
SigmaSS = -L\SigmaSS;
L = expm_fastc(L*dt); %calculations of the Liouvillians


