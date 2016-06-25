% pulse      Pulse definition and excitation profile calculation
%
% [t,y] = pulse(Par)
% [t,y] = pulse(Par,Opt)
%
% [t,y,exprof] = pulse(Par,Opt)
% [t,y,exprof,modulation] = pulse(Par,Opt)
%
% Input: 
%   Par = structure containing the following fields:
%     Par.tp          = pulse length, in us
%     Par.TimeStep    = time step for waveform definition, in us
%                        (default = determined automatically based on
%                         pulse parameters)
%     Par.Flip        = pulse flip angle, in radians (*1) (default: pi)
%     Par.Amplitude   = pulse amplitude, in MHz; ignored if Par.Flip given
%     Par.CenterFreq  = pulse center frequency (default: 0)
%     Par.Phase       = phase for the pulse in radians (default: 0 = +x)
%     Par.Type        = pulse shape name in a string with structure
%                       'AM/FM' (or just 'AM'), where AM refers to the
%                       amplitude modulation function and FM to the
%                       frequency modulation function. The available
%                       options are listed below
%                       If only a single keyword is given, the FM
%                       function is set to 'none'. For a pulse with
%                       constant amplitude, the AM needs to be specified
%                       as 'rectangular'.
%                       Different AM functions can be multiplied by
%                       concatenating the keywords as 'AM1*AM2/FM'.
%                       (default: 'rectangular')
%       Par.*         = value for the pulse parameters defining
%                       the specified pulse shape. The pulse
%                       parameters required for each of the available
%                       modulation functions are listed below.
%       Par.I, Par.Q  = I and Q data describing an arbitrary pulse.
%                       The time axis is reconstructed based on Par.tp
%                       and the length of the I and Q vectors, all
%                       other input parameters (Amplitude, Flip,
%                       CenterFreq, Phase, etc.) are ignored.
%
% Opt = optional structure with the following fields
%       Opt.OverSampleFactor = oversampling factor for the determination
%                              of the time step (default: 10)
%       Opt.Detect           = 'Sz' (default), 'Sy', 'Sx', 'all', define
%                              detection operator for excitation profile
%                              plotting
%       Opt.Offsets          = axis of frequency offsets in MHz for which
%                              to compute the excitation profile
%                              (default approximately ±1.5*BW centered at
%                              Par.CenterFreq, 201 points)
%
% Available pulse modulation functions:
%   - Amplitude modulation: rectangular, gaussian, sinc, quartersin, sech, WURST
%   - Frequency modulation: none, linear, tanh, BWcompensated, uniformQ
%
% The parameters required for the different modulation functions are:
% Amplitude modulation:
% 'rectangular'         - none
% 'gaussian'            - tFWHM     = FWHM in us
%                       Alternatively:
%                       - trunc     = truncation parameter (0 to 1)
% 'sinc'                - zerocross = width between the first zero-
%                                     crossing points in us
% 'sech'                - beta      = dimensionless truncation parameter
%                       - n         = order of the secant function argument
%                                     (default = 1)
% 'WURST'               - nwurst    = WURST n parameter (determining the
%                                     steepness of the amplitude function)
% 'quartersin'          - trise     = rise time in us for quarter sine
%                                     weighting at the pulse edges
%
% Frequency modulation:
% 'linear'              - BW        = frequency sweep width in MHz
%                                     (the frequency is swept over ±BW/2)
%                       - SweepDirection = +1 or -1 (+1 default)
% 'tanh'                - BW        = frequency sweep width in MHz
%                       - beta      = dimensionless truncation parameter
%                       - SweepDirection = +1 or -1 (+1 default)
% 'BWcompensated'(*2)   - BW        = frequency sweep width in MHz
%                                     (the frequency is swept over ±BW/2)
%                       - SweepDirection = +1 or -1 (+1 default)
%                       Parameters required for resonator bandwidth
%                       compensation:
%                       - freqaxis  = frequency axis
%                       - v1        = magnitude response function (ideal or
%                                     experimental)
%                       - Par.mwFreq needs to be defined
% 'uniformQ'            - BW        = frequency sweep width in MHz
%                       The frequency modulation is calculated as the
%                       integral of the squared amplitude modulation
%                       function (for nth order sech/tanh pulses or in
%                       general to obtain offset-independent adiabaticity
%                       pulses, see J. Magn. Reson. 153, 155 (2001).
%
% Output:   t          = time axis for defined waveform in us
%           y          = real and imaginary part of the pulse function
%           If three or four output arguments are requested, the excitation
%           profile is calculated. Additional output fields are:
%           exprof     = structure with frequency offset axis (p.offsets) 
%                        for excitation profile in MHz and excitation
%                        profile (Mi/M0, i = x,y,z, p.Mx, p.My, p.Mz)
%           modulation = structure with amplitude (modulation.A, in MHz),
%                        frequency (modulation.nu, in MHz) and phase
%                        (modulation.phi, in rad) modulation functions
%
% (*1) The conversion from flip angles to amplitudes is performed using the
%      approximations described in:
%      Jeschke, G., Pribitzer, S., Doll, A. Coherence Transfer by Passage
%      Pulses in Electron Paramagnetic Resonance Spectroscopy.
%      J. Phys. Chem. B 119, 13570–13582 (2015).
% (*2) Chirps with variable rate to compensate for the resonator bandwidth.
%      The bandwidth compensation is implemented as described in:
%      Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G., Adiabatic and
%      fast passage ultra-wideband inversion in pulsed EPR
%      J. Magn. Reson. 230, 27-39 (2013).

function varargout = pulse(varargin)

% ----------------------------------------------------------------------- %
% Input argument parsing
% ----------------------------------------------------------------------- %
switch nargin
  case 0
    help(mfilename);
    return
  case 1 % [t,y] = pulse(Par)
    Par = varargin{1};
    Opt = struct;
  case 2 % [t,y] = pulse(Par,Opt)
    Par = varargin{1};
    Opt = varargin{2};
  otherwise
    error('pulse() needs 1 or 2 input arguments.')
end

if ~isstruct(Par)
  error('First input argument must be a structure.');
end
if ~isstruct(Opt)
  error('Second input argument must be a structure.');
end

plotResults = (nargout==0);

if plotResults
  Opt.ExciteProfile = true;
elseif ~isfield(Opt,'ExciteProfile')
  Opt.ExciteProfile = (nargout==3) || (nargout==4);
end

% Set parameters to defaults
%----------------------------------------------------------------------
if ~isfield(Par,'tp')
  error('Pulse length not defined in Par.tp.')
end
if ~isfield(Par,'Flip')
  if ~isfield(Par,'Amplitude')
    Par.Flip = pi;
    Par.Amplitude = [];
  else
    Par.Flip = [];
  end
else
  Par.Amplitude = [];
end

if ~isfield(Par,'CenterFreq')
  Par.CenterFreq = 0; % MHz
end
if ~isfield(Par,'Phase')
  Par.Phase = 0; % rad
end
if ~isfield(Par,'Type')
  Par.Type = 'rectangular';
end

% Options
% ----------------------------------------------------------------------- %
if ~isfield(Opt,'OverSampleFactor')
  Opt.OverSampleFactor = 10;
end
if ~isfield(Opt,'Detect')
  Opt.Detect = 'Sz';
end
if ~isfield(Opt,'nOffsets')
  Opt.nOffsets = 201;
end

% ----------------------------------------------------------------------- %
% Calculate pulse function
% ----------------------------------------------------------------------- %
modulation = struct;
p = struct;

% Check if pulse I and Q data is given
if (isfield(Par,'I') && ~isempty(Par.I)) || ...
   (isfield(Par,'Q') && ~isempty(Par.Q))
  
  if ~isfield(Par,'Type') || isempty(Par.Type)
    Par.Type = 'user-IQ';
  end
  
  if ~isfield(Par,'I')
    Par.I = zeros(size(Par.Q));
  elseif ~isfield(Par,'Q')
    Par.Q = zeros(size(Par.I));
  elseif size(Par.I)~=size(Par.Q);
    error('I and Q input data have different lengths.');
  end
  if ~isvector(Par.I) || ~isvector(Par.Q)
    error('I and Q input data should be vectors.');
  end
  
  t = linspace(0,Par.tp,numel(Par.I));
  nPoints = numel(t);
  y = complex(Par.I,Par.Q);
  
  Par.TimeStep = t(2)-t(1);
  modulation.A = [];
  modulation.dnu = [];
  modulation.phi = [];
  
  if Opt.ExciteProfile && ~isfield(Opt,'Offsets')
    error('For the excitation profile of a user-defined pulse, specify the desired range of frequency offsets Opt.Offsets.');
  end
  
else
  
  % Set pulse shape to rectangular if it is not specified
  if ~isfield(Par,'Type') || isempty(Par.Type)
    Par.Type = 'rectangular';
  end
  
  % Determine pulse shape from input string
  shape = regexp(Par.Type,'/','split');
  switch numel(shape)
    case 1
      AmplitudeModulation = shape;
      FrequencyModulation = 'none';
    case 2
      AmplitudeModulation = regexp(shape{1},'*','split');
      if ~isempty(shape{2})
        FrequencyModulation = shape{2};
      else
        FrequencyModulation = 'none';
      end
    otherwise
      error('Par.Type must contain the pulse shape in the form ''AM/FM''.');
  end
  
  % Check that all the required pulse parameters are given
  for na = 1:numel(AmplitudeModulation)
    switch AmplitudeModulation{na}
      case 'rectangular'
        % no additional parameters needed

      case 'gaussian'
        if (~isfield(Par,'tFWHM') || isempty(Par.tFWHM)) && ...
           (~isfield(Par,'trunc') || isempty(Par.trunc))
          error(['Pulse AM function not sufficiently defined. ',...
            'Specify Par.tFWHM or Par.trunc for the Gaussian envelope.']);
        elseif ~isfield(Par,'tFWHM') && isfield(Par,'trunc') && ...
            ~isempty(Par.trunc)
          % Convert truncation parameter to tFWHM
          Par.tFWHM = sqrt(-(Par.tp^2)/log2(Par.trunc));
          if Par.tFWHM==0
            Par.tFWHM = Par.TimeStep/2;
          end
        end
        
      case 'sinc'
        
        if ~isfield(Par,'zerocross') || isempty(Par.zerocross)
          error(['Pulse AM function not sufficiently defined. ',...
            'Specify Par.zerocross for the sinc envelope.']);
        end
        
      case 'quartersin'
        
        if ~isfield(Par,'trise') || isempty(Par.trise)
          error(['Pulse AM function not sufficiently defined. ',...
            'Specify Par.trise for the quartersine envelope.']);
        end
        
      case 'sech'
        
        if ~isfield(Par,'beta') || isempty(Par.beta)
          error(['Pulse AM function not sufficiently defined. ',...
            'Specify Par.beta parameter (in 1/us) for the sech envelope.']);
        end
        if ~isfield(Par,'n') || isempty(Par.n)
          Par.n = 1;
        end
        
      case 'WURST'
        
        if ~isfield(Par,'nwurst') || isempty(Par.nwurst)
          error(['Pulse AM function not sufficiently defined. ',...
            'Specify Par.nwurst parameter for the WURST envelope.']);
        end
        if numel(Par.nwurst)~=1 || mod(Par.nwurst,1) || Par.nwurst<1
          error('Pulseshape.nwurst must be a nonnegative integer (1,2,...).');
        end
        
      otherwise
        
        error('The amplitude modulation function ''%s'' is not defined.',AmplitudeModulation{na});
        
    end
  end
  
  switch FrequencyModulation
    
    case 'none'
      
    case 'linear'
      
      if ~isfield(Par,'BW') || isempty(Par.BW)
        error(['Pulse FM function not sufficiently defined. ',...
          'Specify Par.BW parameter (in MHz) for the linear chirp.']);
      end
      
    case 'tanh'
      
      if ~isfield(Par,'BW') || isempty(Par.BW)
        error(['Pulse FM function not sufficiently defined. ',...
          'Specify Par.BW parameter (in MHz) for tanh.']);
      end
      if ~isfield(Par,'beta') || isempty(Par.beta)
        error(['Pulse FM function not sufficiently defined. ',...
          'Specify dimensionless Par.beta parameter for tanh.']);
      end
      
    case 'BWcompensated'
      
      if (~isfield(Par,'freqaxis') || isempty(Par.freqaxis)) || ...
          (~isfield(Par,'v1') || isempty(Par.v1))
        error(['Pulse FM function not sufficiently defined. ',...
          'Specify the resonator magnitude response function (Par.freqaxis, Par.v1).']);
      end
      if ~isfield(Par,'mwFreq') || isempty(Par.mwFreq)
        error(['Pulse FM function not sufficiently defined. ',...
          'Par.mwFreq is required to compute resonator bandwidth compensation.']);
      end
      
    case 'uniformQ'
      
      if ~isfield(Par,'BW') || isempty(Par.BW)
        error(['Pulse FM function not sufficiently defined. ',...
          'Specify Par.BW parameter (in MHz).']);
      end
      
    otherwise
      
      error('The frequency modulation function ''%s'' is not defined.',FrequencyModulation);
      
  end
  if ~isfield(Par,'SweepDirection') || isempty(Par.SweepDirection)
    Par.SweepDirection = +1;
  else
    if numel(Par.SweepDirection)~=1 || (Par.SweepDirection~=1 && Par.SweepDirection~=-1)
      error('Par.SweepDirection must be +1 or -1.');
    end
  end
  
  % Set up time axis
  if ~isfield(Par,'TimeStep') || isempty(Par.TimeStep)
    
    % Automatically determine appropriate time step
    %------------------------------------------------------------------ %
    % Determine bandwidth of frequency modulation
    switch FrequencyModulation
      case 'none'
        FM_BW = 0;
      otherwise
        FM_BW = Par.BW;
    end
    
    % Determine bandwidth of amplitude modulation
    AM_BW = 0;
    for na = 1:numel(AmplitudeModulation)
      switch AmplitudeModulation{na}
        case 'rectangular'
          AM_BW = AM_BW + 4/Par.tp;
        case 'gaussian'
          AM_BW = AM_BW + 1/(2*Par.tFWHM);
        case 'sinc'
          AM_BW = AM_BW + 1/Par.zerocross;
        case 'sech'
          % from approximate FWHM of sech function
          AM_BW = AM_BW + Par.beta/(4*Par.tp*asech(0.5));
          %case 'nth order sech'
          % from rise time to 1/2 of maximum value
          %thalf = (Par.tp(n)/2)*(1-(2*asech(0.5)/thisPulse.beta)^(1/thisPulse.n));
          %AM_BW = AM_BW + 1/(4*thalf);
        case 'WURST'
          % rise time from zero to 1/2 of maximum value
          thalf = Par.tp*(1/2 - (1/pi)*asin(2^(-1/Par.nwurst)));
          AM_BW = AM_BW + 1/(4*thalf);
        case 'quartersin'
          AM_BW = AM_BW + 1/(4*Par.trise);
      end
    end
    
    % Calculate maximum frequency offset
    BW = max([FM_BW AM_BW]);
    maxFreq = max(abs(Par.CenterFreq+[-1 1]*BW/2));
    % Use Nyquist theorem to calculate time step, incl. oversampling
    if maxFreq~=0
      Nyquist_dt = 1/(2*maxFreq);
      Par.TimeStep = Nyquist_dt/Opt.OverSampleFactor;
    else
      Par.TimeStep = 0.002; % us
    end
    if Par.TimeStep>Par.tp
      Par.TimeStep = Par.tp;
    end
    Par.TimeStep = Par.tp/round(Par.tp/Par.TimeStep); % last time point = tp
  end
  t = 0:Par.TimeStep:Par.tp;
  ti = t - Par.tp/2;
  nPoints = numel(t);
  
  % ------------------------------------------------------------------- %
  % Amplitude modulation function (modulation.A)
  % ------------------------------------------------------------------- %
  modulation.A = ones(1,nPoints);
  for na = 1:numel(AmplitudeModulation)
    switch AmplitudeModulation{na}
      
      case 'rectangular'
        
        A = ones(1,nPoints);
        
      case 'gaussian'
        
        A = exp(-(4*log(2)*ti.^2)/Par.tFWHM^2);
        
      case 'sinc'
        
        A = sin((2*pi*ti)/Par.zerocross)./((2*pi*ti)/Par.zerocross);
        A(round(Par.tp/(2*Par.TimeStep))+1) = 1;
        A = A/max(A);
        
      case 'quartersin'
        
        % Pulse edges weighted with a quarter period of a sine wave
        A = ones(1,nPoints);
        if Par.trise~=0 && 2*Par.trise<Par.tp
          tpartial = 0:Par.TimeStep:Par.trise;
          npts = numel(tpartial);
          A(1:npts) = sin(tpartial*(pi/(2*Par.trise)));
          A(end-npts+1:end) = A(na,npts:-1:1);
        end
        
      case 'sech'
        
        if Par.n==1 % reduces numerical errors
          A = sech((Par.beta/Par.tp)*ti);
        else
          A = sech(Par.beta*2^(Par.n-1)*(ti/Par.tp).^Par.n);
        end
        
      case 'WURST'
        
        A = 1 - abs(sin(pi*ti/Par.tp)).^Par.nwurst;
        
    end
    modulation.A = modulation.A.*A;
  end
  
  % ------------------------------------------------------------------- %
  % Frequency (modulation.nu) and phase (modulation.phi) modulation functions
  % ------------------------------------------------------------------- %
  switch FrequencyModulation
    
    case 'none'
      
      modulation.nu = zeros(1,nPoints);
      modulation.phi = zeros(1,nPoints);
      
    case 'linear'
      
      k = Par.BW/Par.tp; % frequency sweep (chirp) rate
      modulation.nu = k*ti;
      modulation.phi = 2*pi*Par.SweepDirection*((k/2)*ti.^2);
      
    case 'BWcompensated'
      
      % Variable-rate chirps with resonator bandwidth compensation, as
      % described in:
      %   Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G.,
      %   Adiabatic and fast passage ultra-wideband inversion in
      %   pulsed EPR. J. Magn. Reson. 230, 27–39 (2013).
      %   http://dx.doi.org/10.1016/j.jmr.2013.01.002
      
      % Constant-rate chirp
      k = Par.BW/Par.tp; % frequency sweep rate
      nu_linear = k*ti;
      
      v1_range = interp1(Par.freqaxis,Par.v1,...
        nu_linear+Par.CenterFreq+Par.mwFreq*1e3);
      
      % Frequency dependence of t and time-to-frequency mapping
      c_ = trapz(1./v1_range.^2)/t(end); % const = 2*pi/Qref
      % Qref = reference adiabaticity
      t_f = cumtrapz((1/c_)*v1_range.^-2);
      nu_adapted = interp1(t_f,nu_linear,t,'pchip');
      
      modulation.nu = nu_adapted;
      modulation.phi = 2*pi*Par.SweepDirection*cumtrapz(t,modulation.nu);
      modulation.phi = modulation.phi + abs(min(modulation.phi));  % zero phase offset at pulse center
      
    case 'tanh'
      
      % Determine BWinf parameter from BW and beta parameters
      % (the frequency is swept from -BW/2 to +BW/2)
      Par.BWinf = Par.BW/tanh(Par.beta/2);
      
      modulation.nu = (Par.BWinf/2)*tanh((Par.beta/Par.tp)*ti);
      modulation.phi = (Par.BWinf/2)*(Par.tp/Par.beta)*...
        log(cosh((Par.beta/Par.tp)*ti));
      modulation.phi = 2*pi*Par.SweepDirection*modulation.phi;
      
    case 'uniformQ'
      % The frequency modulation is calculated as the integral of the
      % squared amplitude modulation function (for nth order sech/tanh or
      % in general to obtain offset-independent adiabaticity pulses, see
      %   Garwood, M., DelaBarre, L., J. Magn. Reson. 153, 155-177 (2001)
      %   http://dx.doi.org/10.1006/jmre.2001.2340
      
      modulation.nu = cumtrapz(ti,modulation.A.^2)/trapz(ti,modulation.A.^2); % F2
      modulation.nu = Par.BW*(modulation.nu-1/2);
      modulation.phi = 2*pi*cumtrapz(ti,modulation.nu);
      modulation.phi = modulation.phi + abs(min(modulation.phi)); % zero phase offset at pulse center
      modulation.phi = Par.SweepDirection*modulation.phi;
      
  end
  
  % ------------------------------------------------------------------- %
  % Determine pulse amplitude from flip angle (if only Par.Flip is given)
  % ------------------------------------------------------------------- %
  if (isfield(Par,'Flip') && ~isempty(Par.Flip)) && ...
     (~isfield(Par,'Amplitude') || isempty(Par.Amplitude))
    switch FrequencyModulation
      
      case 'none' % amplitude modulated pulses: beta = integral
        
        Par.Amplitude = Par.Flip/(2*pi*trapz(t,modulation.A));
        
      case {'linear','BWcompensated','tanh','uniformQ'}
        % see
        %    Jeschke et al. (2015) J. Phys. Chem. B, 119, 13570–13582.
        %    http://dx.doi.org/10.1021/acs.jpcb.5b02964
        % Q_crit = (2*pi*v1max)^2/k = minimum adiabaticity on resonance
        
        if Par.Flip>pi
          error('Pulse amplitude calculation from flip angle not applicable for angles larger than pi.')
        end
        Q_crit = (2/pi)*log(2/(1+cos(Par.Flip)));
        if Q_crit>5 % set Q_crit to finite value if it is infinite or very large
          Q_crit = 5;
        end
        
        switch FrequencyModulation
          case {'linear','BWcompensated'}
            sweeprate = Par.BW/Par.tp;
          case 'tanh'
            sweeprate = Par.beta*Par.BWinf/(2*Par.tp);
          case 'uniformQ'
            % Q = w1max^2*A(t)^2/(BW*dnu/dt) see eq. 17 in Garwood, M., DelaBarre, L.,
            % J. Magn. Reson. 153, 155-177 (2001)
            [dummy,ind] = min(abs(ti));
            dnu = diff(2*pi*modulation.nu/(t(2)-t(1)));
            sweeprate = dnu(ind)/(2*pi*(modulation.A(ind))^2);
        end
        
        Par.Amplitude = sqrt(2*pi*Q_crit*sweeprate)/(2*pi);
        
    end
  end
  
  % ------------------------------------------------------------------- %
  % Calculate pulse IQ function
  % ------------------------------------------------------------------- %
  modulation.A = Par.Amplitude*modulation.A;
  totalphase = modulation.phi + 2*pi*Par.CenterFreq*t + Par.Phase;
  y = modulation.A.*exp(1i*totalphase);
  
end

% --------------------------------------------------------------------- %
% Excitation profile calculation
% --------------------------------------------------------------------- %

if Opt.ExciteProfile
  
  % Set up offset axis for excitation profile calculation
  if ~isfield(Opt,'Offsets')
    
    % Estimate bandwidth of the pulse
    if ~strcmp(FrequencyModulation,'none')
      
      BW = Par.BW;
      
    else
      
      % Estimate bandwidth from Fourier transform of amplitude modulation
      % function
      if nextpow2(numel(y))<10
        zf = 2^10;
      else
        zf = 4*2^nextpow2(numel(y));
      end
      yft = abs(fftshift(fft(y,zf)));
      f = fdaxis(Par.TimeStep,zf);
      intg = cumtrapz(yft);
      [dummy,indmax] = min(abs(intg-0.5*max(intg)));
      indbw = find(yft(indmax:end)<0.1*max(yft),1);
      BW = 2*(f(indmax+indbw)-f(indmax));
      
    end
    p.offsets = linspace(-0.75*BW,0.75*BW,Opt.nOffsets) + Par.CenterFreq;
    
  else
    p.offsets = Opt.Offsets;
  end
  nOffsets = numel(p.offsets);
  
  % Spin operators
  Sx = sop(1/2,'x');
  Sy = sop(1/2,'y');
  Sz = sop(1/2,'z');
  
  % Equilibrium density matrix
  Density0 = -Sz;
  
  % Detection operator
  switch Opt.Detect
    case 'Sz', varname = 'Mz';
    case 'Sy', varname = 'My';
    case 'Sx', varname = 'Mx';
    case 'all', 
    otherwise
      error('Opt.Detect has unknown value. Use ''Sx'', ''Sy'', ''Sz'', or ''all''.');
  end
  
  % Pre-allocate result array
  p.Mx = zeros(1,nOffsets);
  p.My = zeros(1,nOffsets);
  p.Mz = zeros(1,nOffsets);
  
  Isignal = real(y);
  Qsignal = imag(y);
  for iOffset = 1:nOffsets
    
    Ham0 = p.offsets(iOffset)*Sz;
    
    % Compute pulse propagator
    if min(y)==max(y) % rectangular pulses
      
      Ham = Isignal(1)*Sx + Qsignal(1)*Sy + Ham0;
      tp = Par.TimeStep*(nPoints-1);
      %UPulse = expm(-2i*pi*Ham*tp);
      M = -2i*pi*tp*Ham; % M = [a b; -b' -a]
      q = sqrt(M(1,1)^2-abs(M(1,2))^2);
      if abs(q)<1e-10
        UPulse = eye(2) + M;
      else
        UPulse = cosh(q)*eye(2) + (sinh(q)/q)*M;
      end
      
    else % general pulses
      
      eye2 = eye(2);
      UPulse = eye2;
      for it = 1:nPoints-1
        
        Ham = Isignal(it)*Sx + Qsignal(it)*Sy + Ham0;
        
        % dU = expm(-2i*pi*Ham*Par.TimeStep);
        % Fast matrix exponential for a traceless, antihermitian 2x2 matrix
        M = -2i*pi*Par.TimeStep*Ham; % M = [a b; -b' -a]
        q = sqrt(M(1)^2-abs(M(3))^2);
        if abs(q)<1e-10
          dU = eye2 + M;
        else
          dU = cosh(q)*eye2 + (sinh(q)/q)*M;
        end
        UPulse = dU*UPulse;
      end
      
    end
    
    % Propagate density matrix
    Density = UPulse*Density0*UPulse';
    
    % Calculate observables
    % (using trace(A*B) = sum(sum(A.*B.')))
    p.Mx(iOffset) = -2*real(sum(sum(Sx.*Density.')));
    p.My(iOffset) = -2*real(sum(sum(Sy.*Density.')));
    p.Mz(iOffset) = -2*real(sum(sum(Sz.*Density.')));
    
  end
  
end

% ----------------------------------------------------------------------- %
% Plotting
% ----------------------------------------------------------------------- %
if plotResults
  clf
  colI = [0 0 1];
  colQ = [1 0 0];
  subplot(2,1,1)
  hold on; box on;
  h = plot(t,modulation.A,t,-modulation.A);
  hI = plot(t,real(y),'Color',colI);
  hQ = plot(t,imag(y),'Color',colQ);
  set(h,'Color',[1 1 1]*0.9);
  Amax = max(modulation.A);
  ylim([-1 1]*Amax*1.1);
  xlabel('t (\mus)')
  ylabel('\nu_1 (MHz)')
  legend([hI hQ],'I','Q')
  axis tight
  title(['Type = ' Par.Type]);

  subplot(2,1,2)
  hold on; box on;
  colBW = [1 1 1]*0.8;
  line([1 1]*Par.CenterFreq,[-1 1],'Color',colBW);
  line([min(p.offsets) max(p.offsets)],[0 0],'Color',colBW);
  if isfield(Par,'BW') && ~isempty(Par.BW) && ~strcmp(FrequencyModulation,'none')
    line([1 1]*(Par.CenterFreq-Par.BW/2),[-1 1],'Color',colBW);
    line([1 1]*(Par.CenterFreq+Par.BW/2),[-1 1],'Color',colBW);
  end
  cc = [0 0 1];
  switch Opt.Detect
    case 'Sz'
      plot(p.offsets,p.Mz,'Color',cc);
      ylabel('{\itM}_z/{\itM}_0')
    case 'Sy'
      plot(p.offsets,p.My,'Color',cc);
      ylabel('{\itM}_y/{\itM}_0')
    case 'Sx'
      plot(p.offsets,p.Mx,'Color',cc);
      ylabel('{\itM}_x/{\itM}_0')
    case 'all'
      plot(p.offsets,p.Mx,...
        p.offsets,p.My,...
        p.offsets,p.Mz);
      ylabel('{\itM}_i/{\itM}_0')
      legend('x','y','z')
  end
  xlabel('frequency (MHz)')
  axis tight
  ylim([-1 1])
end


% ----------------------------------------------------------------------- %
% Output
% ----------------------------------------------------------------------- %
switch nargout
  case 0
    % plotting
  case 2 % [t,y] = pulse(...)
    varargout = {t,y};
  case 3 % [t,y,p] = pulse(...)
    varargout = {t,y,p};
  case 4 % [t,y,p,modulation] = pulse(...)
    varargout = {t,y,p,modulation};
  otherwise
    error('The function pulse() needs 2, 3, or 4 output arguments.')
end
