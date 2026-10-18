% p_lineshape  Convolutional broadening, derivatives and field modulation
%
%   spec = p_lineshape(xAxis,spec,fwhmGL,derivsLeft,mwPhaseLeft,Exp,unitstr,nPadMultiplier,isBroadened)
%
%   Turns an accumulated spectrum into the final detected spectrum.
%
%   Input:
%     xAxis           sweep axis (row vector)
%     spec            spectrum, or one spectrum per row
%     fwhmGL          [Gaussian Lorentzian] convolution FWHM, same units as xAxis
%     derivsLeft      number of derivatives not yet applied
%     mwPhaseLeft     detection phase not yet applied (radians)
%     Exp             experiment structure, processed by p_harmonic; uses
%                     Exp.ModAmp (mT) and Exp.Harmonic for field modulation
%     unitstr         unit of xAxis, for messages
%     nPadMultiplier  range padding for convolution, in units of the larger
%                     FWHM; 0 for no padding
%     isBroadened     true if spec already has line broadening, false if it is
%                     a stick spectrum
%
%   Output:
%     spec            detected spectrum
%
%   Steps:
%   1. Convolution with Lorentzian and Gaussian. The derivatives and the phase
%      are applied with the Lorentzian if it is wide enough, otherwise with the
%      Gaussian. (Both commute with convolution, so the choice does not change
%      the result.) Kernels narrower than 2 increments are skipped.
%   2. Phase via Hilbert transform, if it was not applied in step 1. The
%      dispersion is the Hilbert transform of the absorption, also for sums of
%      lines with different widths (e.g. strains). For a stick spectrum, this
%      is an error.
%   3. Numerical differentiation for any derivatives not applied in step 1.
%      For a stick spectrum, this gives spikes instead of line shapes, so it
%      is an error.
%   4. Field modulation, if Exp.ModAmp>0.

function spec = p_lineshape(xAxis,spec,fwhmGL,derivsLeft,mwPhaseLeft,Exp,unitstr,nPadMultiplier,isBroadened)

dx = xAxis(2) - xAxis(1);
fwhmG = fwhmGL(1);
fwhmL = fwhmGL(2);

% Convolution
%-------------------------------------------------------------------------------
% Kernels narrower than 2 increments have no noticeable effect and are skipped.
% The derivatives and the phase go with the Lorentzian if it is wide enough,
% otherwise with the Gaussian. If neither is wide enough, the phase and the
% derivatives are applied below.
wideL = fwhmL>2*dx;
wideG = fwhmG>2*dx;
if wideL || wideG

  if wideL
    derivL = derivsLeft; phaseL = mwPhaseLeft; derivG = 0; phaseG = 0;
  else
    derivL = 0; phaseL = 0; derivG = derivsLeft; phaseG = mwPhaseLeft;
  end

  % Add padding to left and right of spectral range to reduce convolution artifacts
  nPad = 0;
  if nPadMultiplier>0
    exceedsLowerLimit = any(spec(:,1)~=0);
    exceedsHigherLimit = any(spec(:,end)~=0);
    if exceedsLowerLimit && exceedsHigherLimit
      logmsg(0,'** Spectrum exceeds sweep range. Artifacts at lower and upper limits possible.');
    elseif exceedsLowerLimit
      logmsg(0,'** Spectrum exceeds sweep range. Artifacts at lower limit possible.');
    elseif exceedsHigherLimit
      logmsg(0,'** Spectrum exceeds sweep range. Artifacts at upper limit possible.');
    end
    if exceedsLowerLimit || exceedsHigherLimit
      nPad = round(max(fwhmGL)/dx*nPadMultiplier);
      nPoints = size(spec,2);
      spec = [repmat(spec(:,1),1,nPad) spec repmat(spec(:,end),1,nPad)];
    end
  end

  % Convolution along the sweep dimension (dimension 2 for matrices)
  if isvector(spec), fwhmMask = 1; else, fwhmMask = [0 1]; end

  % Convolution with Lorentzian
  if wideL
    logmsg(1,'  convoluting with Lorentzian, FWHM %g %s, derivative %d',fwhmL,unitstr,derivL);
    spec = convspec(spec,dx,fwhmL*fwhmMask,derivL,0,phaseL);
  end

  % Convolution with Gaussian
  if wideG
    logmsg(1,'  convoluting with Gaussian, FWHM %g %s, derivative %d',fwhmG,unitstr,derivG);
    spec = convspec(spec,dx,fwhmG*fwhmMask,derivG,1,phaseG);
  end

  % Remove padding
  if nPad>0
    spec = spec(:,nPad+(1:nPoints));
  end

  derivsLeft = 0;
  mwPhaseLeft = 0;

end

% Phase via Hilbert transform
%-------------------------------------------------------------------------------
if mwPhaseLeft~=0
  if ~isBroadened
    error(['A dispersion signal (Exp.mwPhase) requires line broadening: a convolutional ' ...
      'line width (Sys.lw) of at least 2 increments (2x%g = %g %s), or strains.'],dx,2*dx,unitstr);
  end
  logmsg(1,'  applying phase %g rad via Hilbert transform',mwPhaseLeft);
  % Linear convolution with the discrete Hilbert kernel 2/(pi*k) for odd k,
  % evaluated via FFT. This keeps the 1/x tails exact over the whole range,
  % unlike the circular FFT-based transform (hilberttrans). Absorption outside
  % the sweep range is missing, which affects Lorentzian lines near the edges.
  nPoints = size(spec,2);
  k = -(nPoints-1):(nPoints-1);
  kernel = zeros(1,2*nPoints-1);
  oddk = mod(k,2)~=0;
  kernel(oddk) = 2./(pi*k(oddk));
  nFFT = 3*nPoints-2;
  specDisp = ifft(fft(spec,nFFT,2).*fft(kernel,nFFT),[],2,'symmetric');
  specDisp = specDisp(:,nPoints:2*nPoints-1);
  spec = cos(mwPhaseLeft)*spec + sin(mwPhaseLeft)*specDisp;
end

% Numerical differentiation
%-------------------------------------------------------------------------------
if derivsLeft>0
  if ~isBroadened
    if ~any(fwhmGL>0)
      error(['Cannot compute a derivative of a stick spectrum without line broadening. ' ...
        'Add a convolutional line width.']);
    else
      error(['Line width (FWHM [%g %g] %s) is smaller than 2 increments (2x%g = %g %s). ' ...
        'Cannot compute derivative. Increase line width or Exp.nPoints.'],fwhmG,fwhmL,unitstr,dx,2*dx,unitstr);
    end
  end
  logmsg(1,'  harmonic %d: using differentiation',derivsLeft);
  for h = 1:derivsLeft
    dspec = diff(spec,[],2)/dx;
    spec = (dspec(:,[1 1:end]) + dspec(:,[1:end end]))/2;
  end
end

% Field modulation
%-------------------------------------------------------------------------------
if Exp.ModAmp>0
  logmsg(1,'  applying field modulation');
  for iSpec = 1:size(spec,1)
    spec(iSpec,:) = fieldmod(xAxis,spec(iSpec,:),Exp.ModAmp,Exp.Harmonic);
  end
end

end

