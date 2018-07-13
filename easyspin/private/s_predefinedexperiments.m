function [Exp] = s_predefinedexperiments(Exp)
% This function creates the Event and Vary Structures

defaultPoints = 512;

P90.tp = 0.016;
P90.Flip = pi/2;

P180.tp = 0.032;
P180.Flip= pi;

switch Exp.Sequence
  case '2pESEEM'
    requiredFields = {'dt','tau'};
    
    checkfields(requiredFields,Exp)
              
    if ~isfield(Exp,'DetWindow')
      if isfield(Exp,'DetSequence')
        Exp = rmfield(Exp.DetSequence);
      end
      Exp.DetWindow = 0;
    end
    
    Exp.Sequence = {P90 Exp.tau P180 Exp.tau};
    
    if ~isfield(Exp,'nPoints')
      Exp.nPoints = defaultPoints;
    end
    
    Exp.Dim1 = {'d1,d2' Exp.dt};
    
  case '3pESEEM'
    
    requiredFields = {'dt','tau','T'};
    
    checkfields(requiredFields,Exp)
        
    if ~isfield(Exp,'DetWindow')
      if isfield(Exp,'DetSequence')
        Exp = rmfield(Exp.DetSequence);
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
    
    requiredFields = {'t1','t2','tau'};
    
    checkfields(requiredFields,Exp)
    
    if ~isfield(Exp,'DetWindow')
      if isfield(Exp,'DetSequence')
        Exp = rmfield(Exp.DetSequence);
      end
      Exp.DetWindow = 0;
    end
    
    Exp.Sequence = {P90 Exp.tau P90 Exp.t1 P180 Exp.t2 P90 Exp.tau};
    
    if ~isfield(Exp,'nPoints')
      Exp.nPoints = [defaultPoints/4 defaultPoints/4];
	  elseif length(Exp.nPoints) == 1
	  Exp.nPoints = [Exp.nPoints Exp.nPoints];
    end
    
    Exp.Dim1 = {'d2' Exp.dt};
    Exp.Dim2 = {'d3' Exp.dt};
    
end


function checkfields(requiredFields, Exp)
    for iField = 1 : numel(requiredFields)
      if ~isfield(Exp,requiredFields{iField})
        msg = ['The field Exp.' requiredFields{iField} ' is required to run ' Exp.Sequence ' with the general method (and may not be allowed to be zero).'];
        error(msg);
      end
    end