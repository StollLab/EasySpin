function SelectedPathways = sf_pathways(Exp)

Display = 0;

%-------------------------------------------------

if (nargin==0)
  Display = 1;
  TestCase = 3;
  switch TestCase
  case 0 % two-pulse ESEEM
    Exp.Inc =     [ +1  +1];
    Exp.t = [200 200]*1e-3;
    Exp.tp = [16 32]*1e-3;
  case 1 % three-pulse ESEEM
    Exp.Inc =     [ 0 +1  0];
    Exp.t = [200 0 200]*1e-3;
    Exp.tp = [16 16 16]*1e-3;
  case 2 % HYSCORE
    Exp.Inc = [ 0 +1  +2  0];
    Exp.t =  [200 0 0 200]*1e-3;
    Exp.tp =  [16 16 16 16]*1e-3;
  case 3
    tau1 = 30; tau2 = 100; T0 = 20;
    Exp.Inc =  [ 0 0 +1 0 0];
    Exp.t =  [tau1 tau1 T0 tau2 tau2]*1e-3;
    Exp.tp = [16 32 16 16 32]*1e-3;
  end
end

nPoints = 256;
dt = 0.08;

nIntervals = numel(Exp.t);

% Generate list of all coherence transfer pathways that are detectable
% (i.e. end in electron coherence order -1).
if nIntervals>1
  N = nIntervals-1;
  [idx{1:N}] = ndgrid(1:4);
  Pathways = fliplr(reshape(cat(ndims(idx{1}),idx{:}),[],N));
  Pathways(:,nIntervals) = 4;
  nPathways = 4^N;
  nEchoPathways = 4^N - 3^N;
  
  CoherenceOrder = [0 0 +1 -1];
  ElectronCoherenceOrder = CoherenceOrder(Pathways);
  
  IncDims = Exp.Inc;
  IntervalTimes = Exp.t;
  
  % Identify all pathways that lead to echoes that either are at the
  % detection point in time or cross it.
  % (Echo is refocused if total time in -1 equals total time in +1)
  for p = 1:nPathways
    
    % Determine number of order changes |delta p|=2
    co = ElectronCoherenceOrder(p,:);
    co(co==0) = [];
    Flips(p) = sum(abs(diff(co))==2);
    
    % Determine echo phase contribution for constant intervals
    t0(p) = ElectronCoherenceOrder(p,:)*IntervalTimes.';
    
    % Determine echo phase contribution for sweep intervals
    if max(abs(IncDims))>0
      
      for d = 1:max(abs(IncDims))
        incdec(p,d) = sum(ElectronCoherenceOrder(p,abs(IncDims)==d));
      end
      
      Refocus(p) = (Flips(p)>0) && (t0(p)==0) && all(incdec(p,:)==0);
      
      Crosses(p) = (Flips(p)>0) && any((t0(p)+nPoints*dt*incdec(p,:))*t0(p)<=0);
      
    else % pulse sequences with echo detection only
      
      Refocus(p) = (Flips(p)>0) && (t0(p)==0);
      
      Crosses(p) = (Flips(p)>0) && (t0(p)==0);
      
    end
    
  end
  
  SelectedPathways = Pathways(Refocus,:);
  
  if (Display)
    % Summary
    fprintf('%d pulses\n%d eCTPs\n%d eCTPs give at least one echo\n',...
      nIntervals,nPathways,nEchoPathways);
    fprintf('%d eCTPs lead to an echo always in the detection window\n',sum(Refocus));
    fprintf('%d crossing echoes\n',sum(Crosses)-sum(Refocus));
    % Symbolic representation of CTPs with +, -, a and b
    Str = 'ab+-';
    fprintf('CTP   #echoes  time    walk dirs\n');
    for iCTP = 1:nPathways
      if Flips(iCTP)>0 && Crosses(iCTP)
        fprintf('%s    %d    %+5.5g     ',Str(Pathways(iCTP,:)),Flips(iCTP),t0(iCTP));
        for iDim=1:max(abs(Exp.Inc))
          fprintf('   %+g',incdec(iCTP,iDim));
        end
        if Refocus(iCTP)
          fprintf('    refocuses');
        elseif Crosses(iCTP)
          fprintf('    crosses');
        end
        fprintf('\n');
      end
    end
    
  end
  
else % FID after single pulse
  
  SelectedPathways = 4;
  
  if (Display)
    % Summary
    fprintf('Single pulse: no echo\n');
  end

end

