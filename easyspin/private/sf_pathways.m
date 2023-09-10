% sf_pathways  Determine refocusing coherence transfer pathways
%
%  pw = sf_pathways(t,inc)
%
% Detemines a list of refocusing pathways given an incrementation scheme
% and initial inter-pulse delays. The analysis assumes zero-length pulses.
%
% Input:
%   t      list of initial inter-pulse delays, in microseconds
%   inc    incrementation scheme, vector of integers
%             0           constant delay
%             +1, +2, ... delay incremented along first, second, etc dimension
%             -1, -2, ... delay decremented along first, second, etc dimension
%
% Output:
%   pw     list of pathways that lead to an echo at the end of the pulse
%          sequence that does not move with delay incrementation
%          1 = alpha, 2 = beta, 3 = +1, 4 = -1

function refocusingPathways = sf_pathways(t,inc)

nPoints = 256;
dt = 0.08;  % µs

Display = nargout==0;

if nargin==0
  Display = true;
  testCase = 3;
  switch testCase
    case 0 % two-pulse ESEEM
      tau0 = 0.2;  % µs
      t = [tau0 tau0];
      inc = [ +1 +1];
    case 1 % three-pulse ESEEM
      tau = 0.2;  % µs
      T0 = 0;
      t = [tau T0 tau];
      inc = [0 +1 0];
    case 2 % HYSCORE
      tau = 0.2;  % µs
      t0 = 0;
      t = [tau t0 t0 tau];
      inc = [0 +1 +2 0];
    case 3
      tau1 = 30;
      tau2 = 100;
      T0 = 20;
      t =  [tau1 tau1 T0 tau2 tau2]*1e-3;
      inc =  [ 0 0 +1 0 0];
  end
end

if numel(inc)~=numel(t)
  error('inc (first input) and t (second input) must have the same number of entries.');
end

nIntervals = numel(t);

if nIntervals>1

  % Generate list of all coherence transfer pathways that are detectable
  % (i.e. end in electron coherence order -1).
  N = nIntervals-1;
  [idx{1:N}] = ndgrid(1:4);
  Pathways = fliplr(reshape(cat(ndims(idx{1}),idx{:}),[],N));
  Pathways(:,nIntervals) = 4;
  nPathways = 4^N;
  nEchoPathways = 4^N - 3^N;
  
  coherenceOrder = [0 0 +1 -1];
  electronCoherenceOrder = coherenceOrder(Pathways);
  
  % Identify all pathways that lead to echoes that either are at the
  % detection point in time or cross it.
  % (Echo is refocused if total time in -1 equals total time in +1)
  for p = 1:nPathways
    
    % Determine number of order changes |delta p|=2
    co = electronCoherenceOrder(p,:);
    co(co==0) = [];
    nFlips(p) = sum(abs(diff(co))==2);
    flipping(p) = nFlips(p)>0;
    
    % Determine echo phase contribution for initial intervals
    t0(p) = electronCoherenceOrder(p,:)*t.';
    
    % Determine echo phase contribution for sweep intervals
    if max(abs(inc))>0

      for d = 1:max(abs(inc))
        incdec(p,d) = sum(electronCoherenceOrder(p,abs(inc)==d));
      end
      refocusing(p) = flipping(p) && (t0(p)==0) && all(incdec(p,:)==0);
      crossing(p) = flipping(p) && any((t0(p)+nPoints*dt*incdec(p,:))*t0(p)<=0);
      
    else % pulse sequences without delay increments
      
      refocusing(p) = flipping(p) && (t0(p)==0);
      crossing(p) = flipping(p) && (t0(p)==0);
      
    end
    
  end
  
  refocusingPathways = Pathways(refocusing,:);
  
  if Display
    % Summary
    fprintf('%d pulses\n%d eCTPs\n%d eCTPs give at least one echo\n',...
      nIntervals,nPathways,nEchoPathways);
    fprintf('%d eCTPs lead to an echo always in the detection window\n',sum(refocusing));
    fprintf('%d crossing echoes\n',sum(crossing)-sum(refocusing));

    % Symbolic representation of CTPs with +, -, a and b
    Str = 'ab+-';
    fprintf('CTP   #echoes  time  type        walk directions\n');
    for iCTP = 1:nPathways
      if nFlips(iCTP)>0 && crossing(iCTP)
        fprintf('%s    %d    %+5.5g',Str(Pathways(iCTP,:)),nFlips(iCTP),t0(iCTP));
        if refocusing(iCTP)
          fprintf('  fixed       ');
        elseif crossing(iCTP)
          fprintf('  crossing    ');
        end
        for iDim = 1:max(abs(inc))
          fprintf('   %+g',incdec(iCTP,iDim));
        end
        fprintf('\n');
      end
    end
    
  end
  
else % FID after single pulse
  
  refocusingPathways = 4;
  
  if Display
    fprintf('Single pulse: no echo\n');
  end

end

end
