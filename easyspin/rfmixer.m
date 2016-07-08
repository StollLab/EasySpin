% rfmixer  Digital up- or downconversion
%
%   [tOut,out] = rfmixer(t,in,LOFreq)     % double-sideband (DSB) mixer
%   [tOut,out] = rfmixer(t,in,LOFreq,'+') % single-sideband (SSB) mixer,
%                                         % upper sideband
%   [tOut,out] = rfmixer(t,in,LOFreq,'-') % single-sideband (SSB) mixer,
%                                         % lower sideband
%
%   [tOut,out] = rfmixer(t,Iin,Qin,LOFreq)               % IQ mixer
%
%   [tOut,Iout,Qout] = rfmixer(t,Iin,Qin,LOFreq)         % IQ modulation
%   [tOut,Iout,Qout] = rfmixer(t,Iin,Qin,LOFreq,'+')
%   [tOut,Iout,Qout] = rfmixer(t,in,LOFreq,'-')          % IQ demodulation
%
%   rfmixer(...)
%   ... = rfmixer(...,Opt)
%
%   Mixes the input signal with the LO frequency. Depending on the number
%   of input and output arguments, the function acts as a double-sideband,
%   single-sideband or IQ mixer or performs IQ modulation or demodulation.
%   Resampling of the input signal before up/downconversion can be
%   performed if the new sampling time step is given.
%
%   If no outputs are requested, the results are plotted.
%
%   Input:
%   - t:             time axis, in ns or us
%   - input signal:  In        = vector of input signal for up/downconversion
%                                (IF or RF) for DSB or SSB mixer or IQ demodulation
%                    OR
%                    Iin,Qin = real and imaginary part of the input signal
%                                for IQ mixer or IQ modulation
%   - LOFreq:        LO frequency, in GHz if t is in ns, or MHz if t is in us
%   - sideband:      absent or '+-' = return both sidebands
%                    '+' = return upper sideband
%                    '-' = return lower sideband
%                    Single-sideband selection ('+' or '-') is only possible
%                    for a real input signal.
%   - options:
%      Opt.dt =    time step for output signal, in ns or us (same as t)
%                  If no time step for resampling is given and the input
%                  time axis is too large, a new time step is computed as
%                  1/(2*Opt.OverSampleFactor*maxfreq).
%      Opt.OverSampleFactor     = factor for determining new time step in
%                                 signal resampling (default = 1.25).
%      Opt.InterpolationMethod  = interpolation method for signal resampling
%                                 (default = 'spline').
%      Opt.BandwidthThreshold   = threshold for input bandwidth determination
%                                 (default = 0.1).
%      Opt.NoiseCutoffThreshold = threshold for selection of part of the input
%                                 signal used to compute amplitude and cos(phase)
%                                 FT overlap (the evaluation of the validity
%                                 of the Hilbert transform is affected by
%                                 the noise level in the signal and baseline).
%      Opt.HilbertThreshold     = threshold for amplitude and phase FT overlap,
%                                 the Hilbert transform corresponds to the
%                                 quadrature signal only if amplitude and phase
%                                 FT do not overlap. If the maximum of the
%                                 product of the normalized amplitude and
%                                 phase FT (extracted from real input signal
%                                 by Hilbert transform) is larger than the
%                                 specified threshold, recovery of the I and
%                                 Q data from the real input signal is not
%                                 possible (default = 0.01).
%
%   Output:
%   - tOut:          time axis for output signal
%   - output signal: out         = up/downconverted signal for DSB, SSB or
%                                  IQ mixer
%                    OR
%                    Iout,Qout = real and imaginary part of the
%                                  up/downconverted signal for IQ modulation
%                                  or demodulation
%
%   Example:
%      t = 0:0.1:200; % ns
%      f = 0.100; % GHz
%      signal_in = cos(2*pi*f*t);
%      LOfreq = 0.5; % GHz
%      rfmixer(t,signal_in,LOfreq,'+');

function varargout = rfmixer(varargin)

if nargin==0
  help(mfilename);
  return
end

% Argument parsing
%-----------------------------------------------------------
t = varargin{1};
switch nargin
  case 3    % rfmixer(t,in,LOFreq) % DSB mixer
    signal = varargin{2};
    LOFreq = varargin{3};
    Opt.sideband = '+-';
    
  case 4
    if ischar(varargin{4}) % rfmixer(t,in,LOFreq,sideband) % SSB
      signal = varargin{2};
      LOFreq = varargin{3};
      Opt.sideband = varargin{4};
    elseif isstruct(varargin{4}) % rfmixer(t,in,LOFreq,Opt) % DSB
      signal = varargin{2};
      LOFreq = varargin{3};
      Opt = varargin{4};
      Opt.sideband = '+-';
    else % [t,out] = rfmixer(t,Iin,Qin,LOFreq) % IQ mixer, IQ modulation
      if ~isequal(size(varargin{2}),size(varargin{3}))
        error('The real and imaginary part of the signal vector must have the same size.')
      end
      signal = varargin{2} + 1i*varargin{3};
      LOFreq = varargin{4};
    end
    
  case 5
    if ischar(varargin{4}) % rfmixer(t,in,LOFreq,sideband,Opt) % SSB
      signal = varargin{2};
      LOFreq = varargin{3};
      Opt = varargin{5};
      Opt.sideband = varargin{4};
    elseif isscalar(varargin{4}) && isstruct(varargin{5}) % rfmixer(t,Iin,Qin,LOFreq,Opt) % IQ
      if ~isequal(size(varargin{2}),size(varargin{3}))
        error('The real and imaginary part of the signal vector must have the same size.')
      end
      signal = varargin{2} + 1i*varargin{3};
      LOFreq = varargin{4};
      Opt = varargin{5};
    elseif ischar(varargin{5}) % rfmixer(t,Iin,Qin,LOfreq,sideband) % IQ
      if ~isequal(size(varargin{2}),size(varargin{3}))
        error('The real and imaginary part of the signal vector must have the same size.')
      end
      signal = varargin{2} + 1i*varargin{3};
      LOFreq = varargin{4};
      Opt.sideband = varargin{5};
    end
    
  case 6
    if ischar(varargin{5}) && isstruct(varargin{6}) % rfmixer(t,Iin,Qin,LOFreq,sideband,Opt) % IQ
      if ~isequal(size(varargin{2}),size(varargin{3}))
        error('The real and imaginary part of the signal vector must have the same size.')
      end
      signal = varargin{2} + 1i*varargin{3};
      LOFreq = varargin{4};
      Opt = varargin{6};
      Opt.sideband = varargin{5};
    end
    
  otherwise
    error('rfmixer() cannot take %d input arguments.',nargin);
    
end

% Input checks
%-----------------------------------------------------------
if ~isfield(Opt,'sideband')
  Opt.sideband = '+-';
end
if ~any(strcmp(Opt.sideband,{'+','-','+-'}))
  error('The sideband specification must be either ''+'', ''-'' or ''+-''.');
end
isSSB = any(strcmp(Opt.sideband,{'+','-'}));

dim = size(signal);
if min(dim)>1 || numel(dim)>2
  error('Input must be a row or column vector.');
end
if numel(t)~=numel(signal)
  error('The signal and its time axis have different lengths.');
end
if size(t)~=size(signal)
  t = reshape(t,size(signal));
end

% Default options
%-----------------------------------------------------------
% Threshold for input bandwidth determination
if ~isfield(Opt,'BandwidthThreshold'), Opt.BandwidthThreshold = 0.1; end 
% Interpolation method for resampling
if ~isfield(Opt,'InterpolationMethod'), Opt.InterpolationMethod = 'spline'; end 
% Factor for oversampling in the determination of the new time step in resampling
if ~isfield(Opt,'OverSampleFactor'), Opt.OverSampleFactor = 1.25; end
% Threshold for part of the input signal used to compute amplitude and
% cos(phase) FT overlap (the evaluation of the validity of the Hilbert
% transform is affected by the noise level in the signal and baseline)
if ~isfield(Opt,'NoiseCutoffThreshold'), Opt.NoiseCutoffThreshold = 0.1; end
% Threshold for amplitude and cos(phase) FT overlap
if ~isfield(Opt,'HilbertThreshold'), Opt.HilbertThreshold = 0.05; end
% Hilbert transform corresponds to quadrature signal only if amplitude
% and cos(phase) FT do not overlap.
% (based on Bedrosian's product theorem; see
%  Boashash, Proc. IEEE 80, 1992, 520, http://dx.doi.org/10.1109/5.135376)

% Signal resampling
%-----------------------------------------------------------
% Determine max frequency of input and output signals
f = fdaxis(t); % GHz or MHz
signalFT = fftshift(fft(signal));
inputband = f(signalFT>Opt.BandwidthThreshold*max(abs(signalFT)));
maxFreqIn = max(inputband);

if LOFreq > maxFreqIn % for upconversion
  maxFreqOut = LOFreq + maxFreqIn;
else % for downconversion
  maxFreqOut = max([LOFreq maxFreqIn]);
end
nyqdt = 1/(2*maxFreqOut); % ns or us, dt for Nyquist criterion

tIn = t;
dtIn = tIn(2) - tIn(1); % input signal time step, in ns or us
if isfield(Opt,'dt')
  if Opt.dt>nyqdt
    fprintf('Warning: Maximum frequency exceeds Nyquist frequency for specified resampling time step.\n')
  end
else
  if dtIn > nyqdt
    Opt.dt = 1/(2*Opt.OverSampleFactor*maxFreqOut); % default resampling time step, in ns or us
  else
    Opt.dt = dtIn;
  end
end
tOut = tIn(1):Opt.dt:tIn(end);
signal_rs = interp1(tIn,signal,tOut,Opt.InterpolationMethod);

% Define I and Q data for real input signals
%-----------------------------------------------------------
singleInput = isreal(signal_rs);
if singleInput
  if isSSB
    
    % Use Hilbert transform to construct quadrature signal
    signal_rs = signal_rs + 1i*imag(hilberttrans(signal_rs));
    
    % Neglect part of the signal amplitude below Opt.NoiseCutoffThreshold
    ind = find(abs(signal_rs)>Opt.NoiseCutoffThreshold*max(abs(signal_rs)));
    signal_nonzero = signal_rs(ind(1):ind(end));

    % Check that input satisfies conditions for Hilbert transform to
    % correspond to the quadrature signal (based on Bedrosian's theorem)
    amplitude_FT = fft(abs(signal_nonzero));
    cosphase_FT = fft(cos(angle(signal_nonzero)));
    amplitude_FT = amplitude_FT/max(amplitude_FT);
    cosphase_FT = cosphase_FT/max(cosphase_FT);
    
    if max(abs(amplitude_FT.*cosphase_FT))>Opt.HilbertThreshold
      error(sprintf(['SSB conversion needs the quadrature signal of the real input signal, which is calculated internally using the Hilbert transform.\n'...
        'In this case, it cannot be recovered. (Use the quadrature signal as input or increase Opt.HilbertThreshold to still perform the calculation.)']));
    end
    
  end
end

% Up/downconversion
%-----------------------------------------------------------
if strcmp(Opt.sideband,'+')
  LOsign = +1;
elseif strcmp(Opt.sideband,'-')
  LOsign = -1;
else
  if LOFreq > maxFreqIn % upconversion
    LOsign = +1;
  else % downconversion
    LOsign = -1;
  end
end
signal_out = signal_rs.*exp(LOsign*2i*pi*LOFreq*tOut);

if singleInput && ~isSSB
  signal_out = real(signal_out);
end

% Output, plotting
%-----------------------------------------------------------
switch nargout
  case 0
    subplot(2,1,1);
    if ~isreal(signal)
      plot(tOut,real(signal_rs),tOut,imag(signal_rs));
      legend('real','imag');
      legend boxoff
    else
      plot(tOut,real(signal_rs));
    end
    title('Input signal');
    xlabel('{\itt}')
    ylim([-1.1 1.1]*max(abs(signal_rs)))
    subplot(2,1,2);
    if ~isreal(signal_out)
      plot(tOut,real(signal_out),tOut,imag(signal_out));
      legend('real','imag');
      legend boxoff
    else
      plot(tOut,real(signal_out));
    end
    title('Output signal');
    xlabel('{\itt}')
    ylim([-1.1 1.1]*max(abs(signal_out)))
    
  case 1
    error('rfmixer() needs to be called with at least two output arguments.\n')
    
  case 2
    varargout{1} = tOut;
    varargout{2} = real(signal_out);
    
  case 3
    varargout{1} = tOut;
    varargout{2} = real(signal_out);
    varargout{3} = imag(signal_out);
   
  otherwise
    error('rfmixer() cannot return %d output arguments.',nargout);
end

return
