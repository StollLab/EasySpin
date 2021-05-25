% validate_exp    Process experimental settings.
%
%   (Could be) Implemented to simplify and maintain consistency in code across programs.
%

function varargout = validate_exp(program,Sys,Exp)

assert(ischar(program), 'Program name must be a string.')

switch program
  case 'pepper'
    % Documented fields and their defaults (mandatory parameters are set to NaN)
    %DefaultExp.mwFreq = NaN; % for field sweeps
    %DefaultExp.Field = NaN; % for frequency sweeps
    DefaultExp.CenterSweep = NaN;
    DefaultExp.Range = NaN;
    DefaultExp.mwCenterSweep = NaN;
    DefaultExp.mwRange = NaN;
    DefaultExp.nPoints = 1024;
    DefaultExp.Temperature = NaN;
    DefaultExp.Harmonic = NaN;
    DefaultExp.Mode = '';
    DefaultExp.Ordering = [];
    DefaultExp.ModAmp = 0;
    DefaultExp.mwPhase = 0;

    DefaultExp.CrystalOrientation = [];
    DefaultExp.CrystalSymmetry = '';
    DefaultExp.MolFrame = [];

    Exp = adddefaults(Exp,DefaultExp);

    % Check microwave frequency and static field
    if ~isfield(Exp,'mwFreq') || isempty(Exp.mwFreq)
      if ~isfield(Exp,'Field')
        error('Please supply either the microwave frequency in Exp.mwFreq (for field sweeps) or the magnetic field in Exp.Field (for frequency sweeps).');
      end
      FieldSweep = false;
    else
      if isfield(Exp,'Field') && ~isempty(Exp.Field)
        error('Give either Exp.mwFreq (for a field sweep) or Exp.Field (for a frequency sweep), but not both.');
      end
      FieldSweep = true;
    end

    if FieldSweep
      if (numel(Exp.mwFreq)~=1) || any(Exp.mwFreq<=0) || ~isreal(Exp.mwFreq)
        error('Uninterpretable microwave frequency in Exp.mwFreq.');
      end
      logmsg(1,'  field sweep, mw frequency %0.8g GHz',Exp.mwFreq);
    else
      if (numel(Exp.Field)~=1) || any(Exp.Field<0) || ~isreal(Exp.Field)
        error('Uninterpretable magnetic field in Exp.Field.');
      end
      logmsg(1,'  frequency sweep, magnetic field %0.8g mT',Exp.Field);
    end

    % Automatic field range determination
    if FieldSweep
      if all(isnan(Exp.CenterSweep)) && all(isnan(Exp.Range))
        if numel(Sys.S)==1 && (Sys.S==1/2) && ~any(Sys.L(:))
          logmsg(1,'  automatic determination of sweep range');
          I = nucspin(Sys.Nucs).';
          if ~isempty(I)
            if Sys.fullA
              Amax = max(abs(Sys.A),[],2);
              Amax = max(reshape(Amax,3,[])).';
            else
              Amax = max(abs(Sys.A),[],2);
            end
          else
            Amax = 0;
          end
          hf = sum(I.*Amax)*1e6; % MHz -> Hz
          if Sys.fullg
            for k = 1:Sys.nElectrons
              g(:,k) = eig(Sys.g((1:3)+(k-1)*3,:));
            end
            g = g(:);
          else
            g = Sys.g(:);
          end
          gmax = max(g);
          gmin = min(g);
          minB = planck*(Exp.mwFreq*1e9 - hf)/bmagn/gmax/1e-3; % mT
          maxB = planck*(Exp.mwFreq*1e9 + hf)/bmagn/gmin/1e-3; % mT
          Center = (maxB+minB)/2; % mT
          Sweep = maxB-minB; % mT
          if Sweep==0, Sweep = 5*max(Sys.lw); end
          if Sweep==0, Sweep = 10; end
          Stretch = 1.25;
          Exp.CenterSweep = [Center, Stretch*Sweep];
        else
          error('Cannot automatically determine field range.\nPlease given either Exp.CenterSweep or Exp.Range.');
        end
      end
    else
      % Automatic range for frequency sweep is done later.
    end

    % Check both CenterSweep and Range, prefer CenterSweep
    if FieldSweep
      if ~isnan(Exp.CenterSweep)
        Exp.Range = Exp.CenterSweep(1) + [-1 1]*Exp.CenterSweep(2)/2;
        Exp.Range = max(Exp.Range,0);
      end
      if isfield(Exp,'Range') && all(~isnan(Exp.Range))
        if any(diff(Exp.Range)<=0) || any(~isfinite(Exp.Range)) || ...
            ~isreal(Exp.Range) || any(Exp.Range<0)
          error('Exp.Range is not valid!');
        end
      end
    else
      if ~isnan(Exp.mwCenterSweep)
        Exp.mwRange = Exp.mwCenterSweep(1) + [-1 1]*Exp.mwCenterSweep(2)/2;
        Exp.mwRange = max(Exp.mwRange,0);
      end
      if isfield(Exp,'mwRange') && all(~isnan(Exp.mwRange))
        if (diff(Exp.mwRange)<=0) || any(~isfinite(Exp.mwRange)) || ...
            ~isreal(Exp.mwRange) || any(Exp.mwRange<0)
          error('Exp.mwRange is not valid!');
        end
      end
    end


    % Number of points
    if any(~isreal(Exp.nPoints)) || numel(Exp.nPoints)>1 || (Exp.nPoints<2)
      error('Problem with Exp.nPoints. Needs to be a number not smaller than 2.')
    end

    if FieldSweep
      logmsg(1,'  frequency %g GHz, field range [%g %g] mT, %d points',...
        Exp.mwFreq,Exp.Range(1),Exp.Range(2),Exp.nPoints);
    else
      if ~SweepAutoRange
        logmsg(1,'  field %g mT, frequency range [%g %g] GHz, %d points',...
          Exp.Field,Exp.mwRange(1),Exp.mwRange(2),Exp.nPoints);
      else
        logmsg(1,'  field %g mT, automatic frequency range, %d points',...
          Exp.Field,Exp.nPoints);
      end
    end

    % Detection harmonic
    autoHarmonic = ~isfield(Exp,'Harmonic') || isempty(Exp.Harmonic) || isnan(Exp.Harmonic);
    noBroadening = (~StrainWidths) && (~ConvolutionBroadening);
    if autoHarmonic
      if FieldSweep && ~noBroadening
        if noBroadening
          Exp.Harmonic = 0;
        else
          Exp.Harmonic = 1;
        end
      else
        Exp.Harmonic = 0;
      end
    end
    if ~any(Exp.Harmonic==[-1,0,1,2])
      error('Exp.Harmonic must be either 0, 1 or 2.');
    end
    if noBroadening && (Exp.Harmonic~=0)
      error('\n  No broadening given. Cannot compute spectrum with Exp.Harmonic=%d.\n',Exp.Harmonic);
    end

    % Modulation amplitude
    if any(Exp.ModAmp<0) || any(isnan(Exp.ModAmp)) || numel(Exp.ModAmp)~=1
      error('Exp.ModAmp must be either a single positive number or zero.');
    end
    if (Exp.ModAmp>0)
      if FieldSweep
        logmsg(1,'  field modulation, amplitude %g mT',Exp.ModAmp);
        if (Exp.Harmonic<1)
          error('With field modulation (Exp.ModAmp), Exp.Harmonic=0 does not work.');
        end
        Exp.ModHarmonic = Exp.Harmonic;
        Exp.ConvHarmonic = 0;
        Exp.DerivHarmonic = 0;
      else
        error('Exp.ModAmp cannot be used with frequency sweeps.');
      end
    else
      Exp.ModHarmonic = 0;
      if ConvolutionBroadening
        Exp.ConvHarmonic = Exp.Harmonic;
        Exp.DerivHarmonic = 0;
      else
        Exp.ConvHarmonic = 0;
        Exp.DerivHarmonic = Exp.Harmonic;
      end
    end

    % Resonator mode
    if ischar(Exp.Mode) && ~isempty(Exp.Mode)
      if strcmp(Exp.Mode,'perpendicular')
      elseif strcmp(Exp.Mode,'parallel')
      else
        error('Exp.Mode must be either ''perpendicular'' or ''parallel''.');
      end
    end
    logmsg(1,'  harmonic %d, %s mode',Exp.Harmonic,Exp.Mode);

    % Powder vs. crystal simulation
    if isfield(Exp,'Orientation') || isfield(Exp,'Orientations')
      error('Exp.Orientation and Exp.Orientations are obsolete (as of EasySpin 5), use Exp.CrystalOrientation instead.');
    end
    PowderSimulation = isempty(Exp.CrystalOrientation);
    Exp.PowderSimulation = PowderSimulation; % for communication with resf*

    % Partial ordering
    if ~isempty(Exp.Ordering)
      if ~PowderSimulation
        error('Partial ordering (Exp.Ordering) can only be used in a powder simulation.');
      end
      if isnumeric(Exp.Ordering) && (numel(Exp.Ordering)==1) && isreal(Exp.Ordering)
        UserSuppliedOrderingFcn = 0;
        logmsg(1,'  partial order (built-in function, lambda = %g)',Exp.Ordering);
      elseif isa(Exp.Ordering,'function_handle')
        UserSuppliedOrderingFcn = 1;
        logmsg(1,'  partial order (user-supplied function)');
      else
        error('Exp.Ordering must be a single number or a function handle.');
      end
      if any(Sys.gStrain) || any(Sys.AStrain) || any(Sys.DStrain) || any(Sys.HStrain)
        error('Exp.Ordering and g/A/D/H strains cannot be used simultaneously.');
      end
    end

    % Temperature and non-equilibrium populations
    NonEquiPops = 0;
    if isfinite(Exp.Temperature)
      if numel(Exp.Temperature)==1
        msg = sprintf('  temperature %g K',Exp.Temperature);
      else
        msg = '  user-specified non-equilibrium populations';
        NonEquiPops = 1;
        if max(Exp.Temperature)==min(Exp.Temperature)
          error('Populations in Exp.Temperature cannot be all equal!');
        end
      end
    else
      msg = '  no temperature';
    end
    
    varargout = {UserSuppliedOrderingFcn, NonEquiPops, msg};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'chili'
    if isfield(Exp,'MOMD')
      error('Exp.MOMD is obsolete. Remove it from your code. See the documentation for details.');
    end

    if ~isfield(Exp,'nPoints'), Exp.nPoints = 1024; end
    if ~isfield(Exp,'Harmonic'), Exp.Harmonic = []; end
    if ~isfield(Exp,'mwPhase'), Exp.mwPhase = 0; end
    if ~isfield(Exp,'Temperature'), Exp.Temperature = NaN; end
    if ~isfield(Exp,'ModAmp'), Exp.ModAmp = 0; end
    if ~isfield(Exp,'Mode'), Exp.Mode = 'perpendicular'; end
    if ~isfield(Exp,'Ordering'), Exp.Ordering = []; end
    if ~isfield(Exp,'CrystalOrientation'), Exp.CrystalOrientation = []; end

    % Number of points
    if any(~isreal(Exp.nPoints)) || numel(Exp.nPoints)>1 || (Exp.nPoints<2)
      error('Problem with Exp.nPoints. Needs to be a number >= 2.')
    end

    % Temperature
    if ~isnan(Exp.Temperature)
      if (numel(Exp.Temperature)~=1) || isinf(Exp.Temperature) || (Exp.Temperature<0)
        error('Problem with Exp.Temperature. If given, Exp.Temperature must be a positive value.')
      end
    end

    % Microwave frequency
    if ~isfield(Exp,'mwFreq')
      if ~isfield(Exp,'Field')
        error('Please supply either the microwave frequency in Exp.mwFreq (for field sweeps) or the magnetic field in Exp.Field (for frequency sweeps).');
      end
      FieldSweep = false;
    else
      if isfield(Exp,'Field')
        error('Give either Exp.mwFreq (for a field sweep) or Exp.Field (for a frequency sweep), but not both.');
      end
      FieldSweep = true;
    end
    if FieldSweep
      if (numel(Exp.mwFreq)~=1) || any(Exp.mwFreq<=0) || ~isreal(Exp.mwFreq)
        error('Uninterpretable microwave frequency in Exp.mwFreq.');
      end
      logmsg(1,'  field sweep, mw frequency %0.8g GHz',Exp.mwFreq);
    else
      if (numel(Exp.Field)~=1) || any(Exp.Field<=0) || ~isreal(Exp.Field)
        error('Uninterpretable magnetic field in Exp.Field.');
      end
      logmsg(1,'  frequency sweep, magnetic field %0.8g mT',Exp.Field);
    end

    % Sweep range (magnetic field, or frequency)
    if FieldSweep
      if isfield(Exp,'CenterSweep')
        if isfield(Exp,'Range')
          logmsg(0,'Using Exp.CenterSweep and ignoring Exp.Range.');
        end
      else
        if isfield(Exp,'Range')
          Exp.CenterSweep = [mean(Exp.Range) diff(Exp.Range)];
        else
          if (Sys.nElectrons==1) && (Sys.S==1/2)
            logmsg(1,'  automatic determination of sweep range');
            Stretch = 1.25;
            I = nucspin(Sys.Nucs).';
            if numel(I)>0
              Amax = max(abs(Sys.A),[],2);
              hf = sum(I.*Amax)*1e6; % MHz -> Hz
            else
              hf = 0;
            end
            gmax = max(Sys.g(:));
            gmin = min(Sys.g(:));
            if FieldSweep
              minB = planck*(Exp.mwFreq*1e9 - hf)/bmagn/gmax/1e-3;
              maxB = planck*(Exp.mwFreq*1e9 + hf)/bmagn/gmin/1e-3;
              Exp.CenterSweep = [(maxB+minB)/2, Stretch*max(maxB-minB,5)];
            else
              minE = bmagn*Exp.Field*1e-3*gmin/planck - hf; % Hz
              maxE = bmagn*Exp.Field*1e-3*gmax/planck + hf; % Hz
              Exp.CenterSweep = [(maxE+minE)/2, Stretch*max(maxE-minE,10e6)]/1e9; % GHz
            end
          else
            error('Cannot automatically determine sweep range for this spin system.');
          end
        end
      end
    else
      if isfield(Exp,'mwCenterSweep')
        if isfield(Exp,'mwRange')
          logmsg(0,'Using Exp.mwCenterSweep and ignoring Exp.mwRange.');
        end
      else
        if isfield(Exp,'mwRange')
          Exp.mwCenterSweep = [mean(Exp.mwRange) diff(Exp.mwRange)];
        else
          error('Either Exp.mwRange or Exp.mwCenterSweep need to be given.');
        end
      end
    end

    if FieldSweep
      CenterField = Exp.CenterSweep(1);
      Sweep = Exp.CenterSweep(2);
      Exp.Range = Exp.CenterSweep(1) + [-1 1]/2*Sweep;
      if any(Exp.Range<0) || diff(Exp.Range)<=0
        error('Invalid sweep range! Check Exp.CenterSweep or Exp.Range.');
      end
    else
      CenterFreq = Exp.mwCenterSweep(1);
      Sweep = Exp.mwCenterSweep(2);
      Exp.mwRange = Exp.mwCenterSweep(1) + [-1 1]/2*Sweep;
      CenterField = Exp.Field;
      if any(Exp.mwRange<0) || diff(Exp.mwRange)<=0
        error('Invalid sweep range! Check Exp.mwCenterSweep or Exp.mwRange.');
      end
    end

    if FieldSweep
      logmsg(1,'  field range (mT): min %g, max %g, center %g, width %g',...
        Exp.Range(1),Exp.Range(2),CenterField,Sweep);
    else
      logmsg(1,'  frequency range (GHz): min %g, max %g, center %g, width %g',...
        Exp.mwRange(1),Exp.mwRange(2),CenterFreq,Sweep);
    end

    % Detection harmonic
    if ~isfield(Exp,'Harmonic') || isempty(Exp.Harmonic) || isnan(Exp.Harmonic)
      if FieldSweep
        Exp.Harmonic = 1;
      else
        Exp.Harmonic = 0;
      end
    end
    if ~any(Exp.Harmonic==[-1,0,1,2])
      error('Exp.Harmonic must be 0, 1 or 2.');
    end

    % Modulation amplitude
    if any(Exp.ModAmp<0) || any(isnan(Exp.ModAmp)) || numel(Exp.ModAmp)~=1
      error('Exp.ModAmp must be either a single positive number or zero.');
    end
    if (Exp.ModAmp>0)
      if FieldSweep
        logmsg(1,'  field modulation, amplitude %g mT',Exp.ModAmp);
        if (Exp.Harmonic<1)
          error('With field modulation (Exp.ModAmp), Exp.Harmonic=0 does not work.');
        end
        Exp.ModHarmonic = Exp.Harmonic;
        Exp.ConvHarmonic = 0;
        Exp.DerivHarmonic = 0;
      else
        error('Exp.ModAmp cannot be used with frequency sweeps.');
      end
    else
      Exp.ModHarmonic = 0;
      if ConvolutionBroadening
        Exp.ConvHarmonic = Exp.Harmonic;
        Exp.DerivHarmonic = 0;
      else
        Exp.ConvHarmonic = 0;
        Exp.DerivHarmonic = Exp.Harmonic;
      end
    end

    % Resonator mode
    switch Exp.Mode
      case 'perpendicular', ParallelMode = false;
      case 'parallel', ParallelMode = true;
      otherwise, error('Exp.Mode must be either ''perpendicular'' or ''parallel''.');
    end
    logmsg(1,'  harmonic %d, %s mode',Exp.Harmonic,Exp.Mode);

    % Complain if fields only valid in pepper() are given
    if isfield(Exp,'Orientations')
      warning('Exp.Orientations is obsolete. Use Exp.CrystalOrientations instead.');
    end
    if isfield(Exp,'CrystalSymmetry')
      warning('Exp.CrystalSymmetry is not used by chili.');
    end

    % Partial ordering
    if ~isempty(Exp.Ordering)
      %if ~PowderSimulation
      %  error('Partial ordering (Exp.Ordering) can only be used in a powder simulation.');
      %end
      if isnumeric(Exp.Ordering) && (numel(Exp.Ordering)==1) && isreal(Exp.Ordering)
        UserSuppliedOrderingFcn = false;
        logmsg(1,'  partial order (built-in function, coefficient = %g)',Exp.Ordering);
      elseif isa(Exp.Ordering,'function_handle')
        UserSuppliedOrderingFcn = true;
        logmsg(1,'  partial order (user-supplied function)');
      else
        error('Exp.Ordering must be a single number or a function handle.');
      end
    end

    % Determine whether to do a powder simulation
    if ~usePotential
      if isempty(Exp.Ordering) || all(Exp.Ordering==0)
        logmsg(1,'  No orientational potential given, skipping powder simulation.');
        PowderSimulation = false;
      else
      logmsg(1,'  Orientational potential given, doing powder simulation.');
        PowderSimulation = true;
      end    
    else
      if ~isempty(Exp.CrystalOrientation)
        logmsg(1,'  Orientational potential given, doing single-crystal simulation.');
        PowderSimulation = false;
      else
        logmsg(1,'  Orientational potential given, doing powder simulation.');
        PowderSimulation = true;
      end
    end
    
    varargout = {Exp, CenterField, ParallelMode, UserSuppliedOrderingFcn, PowderSimulation};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'cardamom'

    if ~isfield(Exp,'nPoints'), Exp.nPoints = 1024; end
    if ~isfield(Exp,'Harmonic'), Exp.Harmonic = []; end
    if ~isfield(Exp,'mwPhase'), Exp.mwPhase = 0; end
    if ~isfield(Exp,'Temperature'), Exp.Temperature = NaN; end
    if ~isfield(Exp,'ModAmp'), Exp.ModAmp = 0; end
    if ~isfield(Exp,'Mode'), Exp.Mode = 'perpendicular'; end
    if ~isfield(Exp,'Ordering'), Exp.Ordering = []; end
    if ~isfield(Exp,'CrystalOrientation'), Exp.CrystalOrientation = []; end

    % Number of points
    if any(~isreal(Exp.nPoints)) || numel(Exp.nPoints)>1 || (Exp.nPoints<2)
      error('Problem with Exp.nPoints. Needs to be a number >= 2.')
    end

%     % Temperature  TODO implement in cardamom
%     if ~isnan(Exp.Temperature)
%       if (numel(Exp.Temperature)~=1) || isinf(Exp.Temperature) || (Exp.Temperature<0)
%         error('Problem with Exp.Temperature. If given, Exp.Temperature must be a positive value.')
%       end
%     end

    % Microwave frequency
    if ~isfield(Exp,'mwFreq')
      if ~isfield(Exp,'Field')
        error('Please supply either the microwave frequency in Exp.mwFreq (for field sweeps) or the magnetic field in Exp.Field (for frequency sweeps).');
      end
      FieldSweep = false;
    else
      if isfield(Exp,'Field')
        error('Give either Exp.mwFreq (for a field sweep) or Exp.Field (for a frequency sweep), but not both.');
      end
      FieldSweep = true;
    end
    if FieldSweep
      if (numel(Exp.mwFreq)~=1) || any(Exp.mwFreq<=0) || ~isreal(Exp.mwFreq)
        error('Uninterpretable microwave frequency in Exp.mwFreq.');
      end
      logmsg(1,'  field sweep, mw frequency %0.8g GHz',Exp.mwFreq);
    else
      if (numel(Exp.Field)~=1) || any(Exp.Field<=0) || ~isreal(Exp.Field)
        error('Uninterpretable magnetic field in Exp.Field.');
      end
      logmsg(1,'  frequency sweep, magnetic field %0.8g mT',Exp.Field);
    end

    % Sweep range (magnetic field, or frequency)
    if FieldSweep
      if isfield(Exp,'CenterSweep')
        if isfield(Exp,'Range')
          logmsg(0,'Using Exp.CenterSweep and ignoring Exp.Range.');
        end
      else
        if isfield(Exp,'Range')
          Exp.CenterSweep = [mean(Exp.Range) diff(Exp.Range)];
        else
          if (Sys.nElectrons==1) && (Sys.S==1/2)
            logmsg(1,'  automatic determination of sweep range');
            Stretch = 1.25;
            I = nucspin(Sys.Nucs).';
            if numel(I)>0
              Amax = max(abs(Sys.A),[],2);
              hf = sum(I.*Amax)*1e6; % MHz -> Hz
            else
              hf = 0;
            end
            gmax = max(Sys.g(:));
            gmin = min(Sys.g(:));
            if FieldSweep
              minB = planck*(Exp.mwFreq*1e9 - hf)/bmagn/gmax/1e-3;
              maxB = planck*(Exp.mwFreq*1e9 + hf)/bmagn/gmin/1e-3;
              Exp.CenterSweep = [(maxB+minB)/2, Stretch*max(maxB-minB,5)];
            else
              minE = bmagn*Exp.Field*1e-3*gmin/planck - hf; % Hz
              maxE = bmagn*Exp.Field*1e-3*gmax/planck + hf; % Hz
              Exp.CenterSweep = [(maxE+minE)/2, Stretch*max(maxE-minE,10e6)]/1e9; % GHz
            end
          else
            error('Cannot automatically determine sweep range for this spin system.');
          end
        end
      end
    else
      if isfield(Exp,'mwCenterSweep')   %TODO implement in cardamom
        if isfield(Exp,'mwRange')
          logmsg(0,'Using Exp.mwCenterSweep and ignoring Exp.mwRange.');
        end
      else
        if isfield(Exp,'mwRange')
          Exp.mwCenterSweep = [mean(Exp.mwRange) diff(Exp.mwRange)];
        else
          error('Either Exp.mwRange or Exp.mwCenterSweep need to be given.');
        end
      end
    end

    if FieldSweep
      CenterFreq = [];
      CenterField = Exp.CenterSweep(1);
      Sweep = Exp.CenterSweep(2);
      Exp.Range = Exp.CenterSweep(1) + [-1 1]/2*Sweep;
      if any(Exp.Range<0) || diff(Exp.Range)<=0
        error('Invalid sweep range! Check Exp.CenterSweep or Exp.Range.');
      end
    else
      CenterFreq = Exp.mwCenterSweep(1);
      Sweep = Exp.mwCenterSweep(2);
      Exp.mwRange = Exp.mwCenterSweep(1) + [-1 1]/2*Sweep;
      CenterField = Exp.Field;
      if any(Exp.mwRange<0) || diff(Exp.mwRange)<=0
        error('Invalid sweep range! Check Exp.mwCenterSweep or Exp.mwRange.');
      end
    end

    if FieldSweep
      logmsg(1,'  field range (mT): min %g, max %g, center %g, width %g',...
        Exp.Range(1),Exp.Range(2),CenterField,Sweep);
    else
      logmsg(1,'  frequency range (GHz): min %g, max %g, center %g, width %g',...
        Exp.mwRange(1),Exp.mwRange(2),CenterFreq,Sweep);
    end

%     % Detection harmonic  TODO implement in cardamom
%     if ~isfield(Exp,'Harmonic') || isempty(Exp.Harmonic) || isnan(Exp.Harmonic)
%       if FieldSweep
%         Exp.Harmonic = 1;
%       else
%         Exp.Harmonic = 0;
%       end
%     end
%     if ~any(Exp.Harmonic==[-1,0,1,2])
%       error('Exp.Harmonic must be 0, 1 or 2.');
%     end

%     % Modulation amplitude  TODO implement in cardamom
%     if any(Exp.ModAmp<0) || any(isnan(Exp.ModAmp)) || numel(Exp.ModAmp)~=1
%       error('Exp.ModAmp must be either a single positive number or zero.');
%     end
%     if (Exp.ModAmp>0)
%       if FieldSweep
%         logmsg(1,'  field modulation, amplitude %g mT',Exp.ModAmp);
%         if (Exp.Harmonic<1)
%           error('With field modulation (Exp.ModAmp), Exp.Harmonic=0 does not work.');
%         end
%         Exp.ModHarmonic = Exp.Harmonic;
%         Exp.ConvHarmonic = 0;
%         Exp.DerivHarmonic = 0;
%       else
%         error('Exp.ModAmp cannot be used with frequency sweeps.');
%       end
%     else
%       Exp.ModHarmonic = 0;
%       if ConvolutionBroadening  % TODO perhaps put declaration in a better spot in chili
%         Exp.ConvHarmonic = Exp.Harmonic;
%         Exp.DerivHarmonic = 0;
%       else
%         Exp.ConvHarmonic = 0;
%         Exp.DerivHarmonic = Exp.Harmonic;
%       end
%     end

%     % Resonator mode  TODO implement in cardamom
%     switch Exp.Mode
%       case 'perpendicular', ParallelMode = false;
%       case 'parallel', ParallelMode = true;
%       otherwise, error('Exp.Mode must be either ''perpendicular'' or ''parallel''.');
%     end
%     logmsg(1,'  harmonic %d, %s mode',Exp.Harmonic,Exp.Mode);
% 
%     % Complain if fields only valid in pepper() are given
%     if isfield(Exp,'Orientations')
%       warning('Exp.Orientations is obsolete. Use Exp.CrystalOrientations instead.');
%     end
%     if isfield(Exp,'CrystalSymmetry')
%       warning('Exp.CrystalSymmetry is not used by chili.');
%     end

%     % Partial ordering  TODO implement in cardamom
%     if ~isempty(Exp.Ordering)
%       %if ~PowderSimulation
%       %  error('Partial ordering (Exp.Ordering) can only be used in a powder simulation.');
%       %end
%       if isnumeric(Exp.Ordering) && (numel(Exp.Ordering)==1) && isreal(Exp.Ordering)
%         UserSuppliedOrderingFcn = false;
%         logmsg(1,'  partial order (built-in function, coefficient = %g)',Exp.Ordering);
%       elseif isa(Exp.Ordering,'function_handle')
%         UserSuppliedOrderingFcn = true;
%         logmsg(1,'  partial order (user-supplied function)');
%       else
%         error('Exp.Ordering must be a single number or a function handle.');
%       end
%     end

%     % Determine whether to do a powder simulation  TODO implement in cardamom
%     if ~usePotential
%       if isempty(Exp.Ordering) || all(Exp.Ordering==0)
%         logmsg(1,'  No orientational potential given, skipping powder simulation.');
%         PowderSimulation = false;
%       else
%       logmsg(1,'  Orientational potential given, doing powder simulation.');
%         PowderSimulation = true;
%       end    
%     else
%       if ~isempty(Exp.CrystalOrientation)
%         logmsg(1,'  Orientational potential given, doing single-crystal simulation.');
%         PowderSimulation = false;
%       else
%         logmsg(1,'  Orientational potential given, doing powder simulation.');
%         PowderSimulation = true;
%       end
%     end
    
    varargout = {Exp, FieldSweep, CenterField, CenterFreq, Sweep};%, ParallelMode, UserSuppliedOrderingFcn, PowderSimulation};

  otherwise
    error('Program not recognized.')
    
end

end