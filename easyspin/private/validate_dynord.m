% check_dynord    Process dynamics and orientational potential.
%
%   (Could be) Implemented to simplify and maintain consistency in code across programs.
%

function varargout = validate_dynord(program,Sys,FieldSweep,isDiffSim)

assert(ischar(program), 'Program name must be a string.')

if nargin<3
  FieldSweep = false;
  isDiffSim = true;
elseif nargin<4
  isDiffSim = true;
end

switch program
  case 'chili'
    if isfield(Sys,'psi')
      error('Sys.psi is obsolete. Remove it from your code. See the documentation for details.');
    end

    if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end
    if ~isfield(Sys,'Exchange'), Sys.Exchange = 0; end
    if ~isfield(Sys,'Potential'), Sys.Potential = []; end

    if isfield(Sys,'tcorr'), Dynamics.tcorr = Sys.tcorr; end
    if isfield(Sys,'Diff'), Dynamics.Diff = Sys.Diff; end
    if isfield(Sys,'logtcorr'), Dynamics.logtcorr = Sys.logtcorr; end
    if isfield(Sys,'logDiff'), Dynamics.logDiff = Sys.logDiff; end
    
    if isfield(Sys,'lwpp'), Dynamics.lwpp = Sys.lwpp; end
    if isfield(Sys,'lw'), Dynamics.lw = Sys.lw; end

    Dynamics.Exchange = Sys.Exchange;
    usePotential = ~isempty(Sys.Potential) && ~all(Sys.Potential(:,4)==0);
    
    [Dynamics,err] = processdynamics(Dynamics,FieldSweep);
    error(err);
    
    varargout = {Dynamics,Potential,usePotential};
    
  case 'cardamom'
    
%     if isfield(Sys,'psi')
%       error('Sys.psi is obsolete. Remove it from your code. See the documentation for details.');
%     end
% 
%     if ~isfield(Sys,'DiffFrame'), Sys.DiffFrame = [0 0 0]; end  % TODO implement in cardamom
%     if ~isfield(Sys,'Exchange'), Sys.Exchange = 0; end

    if ~isfield(Sys,'DiffGlobal'), Sys.DiffGlobal = []; end

    if isDiffSim
      if isfield(Sys,'tcorr'), Dynamics.tcorr = Sys.tcorr; end  % TODO process and feed to stochtraj?
      if isfield(Sys,'Diff'), Dynamics.Diff = Sys.Diff; end
      if isfield(Sys,'logtcorr'), Dynamics.logtcorr = Sys.logtcorr; end
      if isfield(Sys,'logDiff'), Dynamics.logDiff = Sys.logDiff; end
    end
    
    if isfield(Sys,'DiffGlobal'), Dynamics.DiffGlobal = Sys.DiffGlobal; end
    
    if isfield(Sys,'lwpp'), Dynamics.lwpp = Sys.lwpp; end
    if isfield(Sys,'lw'), Dynamics.lw = Sys.lw; end
    
    [Dynamics,err] = processdynamics(Dynamics,FieldSweep,isDiffSim);
    error(err);
    
    varargout = {Dynamics};%,Potential,usePotential}

  case 'stochtraj_diffusion'
    if isfield(Sys,'Potential')
      if isfield(Sys.Potential,'lambda') && isfield(Sys.Potential,'LMK')
        if ~ismatrix(Sys.Potential.LMK) || size(Sys.Potential.LMK,2)~=3
          error('LMK must be an array of shape Nx3.')
        end
        if ~ismatrix(Sys.Potential.lambda) || size(Sys.Potential.lambda,2)~=2
          error('lambda must be an array of shape Nx2.')
        end
        % Enforce indexing convention
        for j=1:size(Sys.Potential.LMK,1)
          L = Sys.Potential.LMK(j,1);
          M = Sys.Potential.LMK(j,2);
          K = Sys.Potential.LMK(j,3);
          assert(L>0,'For all sets of indices LMK, it is required that L>0.')
          if K==0
            assert((0<=M)&&(M<=L),'For all sets of indices LMK, if K=0, then it is required that 0<=M<=L.')
          else
            assert((0<K)&&(K<=L)&&abs(M)<=L,'For all sets of indices LMK, if K~=0, then it is required that 0<K<=L and |M|<=L.')
          end
        end
        Sim.lambda = Sys.Potential.lambda;
        Sim.LMK = Sys.Potential.LMK;
      elseif ~isfield(Sys.Potential,'lambda') && ~isfield(Sys.Potential,'LMK')
        % if no orientational potential coefficient is given, initialize empty arrays
        Sim.lambda = [];
        Sim.LMK = [];
      else
        error('Both potential coefficients and LMK are required for an orientational potential.')
      end
    else
      Sim.lambda = [];
      Sim.LMK = [];
    end

    % parse the dynamics parameter input using private function
    if isfield(Sys,'tcorr'), Dynamics.tcorr = Sys.tcorr; end
    if isfield(Sys,'Diff'), Dynamics.Diff = Sys.Diff; end
    if isfield(Sys,'logtcorr'), Dynamics.logtcorr = Sys.logtcorr; end
    if isfield(Sys,'logDiff'), Dynamics.logDiff = Sys.logDiff; end
    
    % FieldSweep not implemented for stochtraj yet
    [Dynamics, err] = processdynamics(Dynamics,[],isDiffSim);
    error(err);

    varargout = {Dynamics,Sim};

  otherwise
    error('Program not recognized.')
    
end

end


% Helper function
% -------------------------------------------------------------------------
function [Dyn,err] = processdynamics(D,FieldSweep,isDiffSim)

Dyn = D;
err = '';

if isDiffSim
  % diffusion tensor, correlation time
  %------------------------------------------------------------------------
  % convert everything (tcorr, logcorr, logDiff) to Diff
  if isfield(Dyn,'Diff')
    % Diff given
  elseif isfield(Dyn,'logDiff')
    Dyn.Diff = 10.^Dyn.logDiff;
  elseif isfield(Dyn,'tcorr')
    Dyn.Diff = 1/6./Dyn.tcorr;
  elseif isfield(Dyn,'logtcorr')
    if Dyn.logtcorr>=0, error('Sys.logtcorr must be negative.'); end
    Dyn.Diff = 1/6./10.^Dyn.logtcorr;
  else
    err = sprintf('You must specify a rotational correlation time or a diffusion tensor\n(Sys.tcorr, Sys.logtcorr, Sys.Diff or Sys.logDiff).');
    return
  end

  if any(Dyn.Diff<0)
    error('Negative diffusion rate or correlation times are not possible.');
  elseif any(Dyn.Diff>1e12)
    fprintf('Diffusion rate very fast. Simulation might not converge.\n');
  elseif any(Dyn.Diff<1e3)
    fprintf('Diffusion rate very slow. Simulation might not converge.\n');
  end

  % expand to rhombic tensor
  switch numel(Dyn.Diff)
    case 1, Dyn.Diff = Dyn.Diff([1 1 1]);
    case 2, Dyn.Diff = Dyn.Diff([1 1 2]);
    case 3, % Diff already rhombic
    otherwise
      err = 'Sys.Diff must have 1, 2 or 3 elements (isotropic, axial, rhombic).';
      return
  end

end

if isfield(Dyn,'lw')
  if numel(Dyn.lw)>1
    if FieldSweep
      LorentzFWHM = Dyn.lw(2)*28 * 1e6; % mT -> MHz -> Hz
    else
      LorentzFWHM = Dyn.lw(2)*1e6; % MHz -> Hz
    end
  else
    LorentzFWHM = 0;
  end
  if LorentzFWHM~=0
    % Lorentzian T2 from FWHM in freq domain 1/T2 = pi*FWHM
    Dyn.T2 = 1/LorentzFWHM/pi;
  else
    Dyn.T2 = Inf;
  end
end

% Heisenberg exchange
%------------------------------------------------------------------
if ~isfield(Dyn,'Exchange'), Dyn.Exchange = 0; end
Dyn.Exchange = Dyn.Exchange*2*pi*1e6; % MHz -> angular frequency

end