function [Exp] = s_predefinedexperiments(Exp)
% This function creates the Event and Vary Structures

defaultPoints = 512;

% Predefined experiments use pulses along y, as in saffron's fast mode
P90.tp = 0.016;
P90.Flip = pi/2;
P90.Phase = pi/2;

P180.tp = 0.032;
P180.Flip = pi;
P180.Phase = pi/2;

if ~ismember(Exp.Sequence,{'2pESEEM','3pESEEM','HYSCORE'})
  error('Exp.Sequence = ''%s'' is not supported with the thyme method.',Exp.Sequence);
end

if ~isfield(Exp,'dt')
  error('Exp.dt is missing.');
end

% Initial values of incremented delays default to zero, as in saffron's fast mode
switch Exp.Sequence
  case '2pESEEM'
    if ~isfield(Exp,'tau'), Exp.tau = 0; end
  case '3pESEEM'
    if ~isfield(Exp,'T'), Exp.T = 0; end
  case 'HYSCORE'
    if ~isfield(Exp,'t1'), Exp.t1 = 0; end
    if ~isfield(Exp,'t2'), Exp.t2 = 0; end
end
if ismember(Exp.Sequence,{'3pESEEM','HYSCORE'})
  if ~isfield(Exp,'tau')
    error('Exp.tau is missing.');
  end
  if Exp.tau<=0
    error('Exp.tau must be larger than 0.');
  end
end

switch Exp.Sequence
  case '2pESEEM'

    if ~isfield(Exp,'DetWindow')
      if isfield(Exp,'DetSequence')
        Exp = rmfield(Exp,'DetSequence');
      end
      Exp.DetWindow = 0;
    end
    
    Exp.Sequence = {P90 Exp.tau P180 Exp.tau};
    
    if ~isfield(Exp,'nPoints')
      Exp.nPoints = defaultPoints;
    end
    
    Exp.Dim1 = {'d1,d2' Exp.dt};
    
  case '3pESEEM'
    
    if ~isfield(Exp,'DetWindow')
      if isfield(Exp,'DetSequence')
        Exp = rmfield(Exp,'DetSequence');
      end
      Exp.DetWindow = 0;
    end
    
    Exp.Sequence = {P90 Exp.tau P90 Exp.T P90 Exp.tau};
    
    if ~isfield(Exp,'nPoints')
      Exp.nPoints = defaultPoints;
    end
    
    Exp.Dim1 = {'d2' Exp.dt};
    
    Exp.PhaseCycle{1}=[0, 1; pi, -1];
    Exp.PhaseCycle{3}=[0, 1; pi, -1];

  case 'HYSCORE'
    
    if ~isfield(Exp,'DetWindow')
      if isfield(Exp,'DetSequence')
        Exp = rmfield(Exp,'DetSequence');
      end
      Exp.DetWindow = 0;
    end
    
    Exp.Sequence = {P90 Exp.tau P90 Exp.t1 P180 Exp.t2 P90 Exp.tau};
    
    if ~isfield(Exp,'nPoints')
      Exp.nPoints = [defaultPoints/4 defaultPoints/4];
	  elseif length(Exp.nPoints) == 1
	  Exp.nPoints = [Exp.nPoints Exp.nPoints];
    end
    
    % one increment for both dimensions, or one per dimension
    if numel(Exp.dt)>2
      error('Exp.dt needs either 1 or 2 elements, one per dimension. You gave %d.',numel(Exp.dt));
    end
    dt = Exp.dt.*[1 1];
    Exp.Dim1 = {'d2' dt(1)};
    Exp.Dim2 = {'d3' dt(2)};

end
