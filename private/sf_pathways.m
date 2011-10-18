function SelectedCTP = sf_pathways(Exp)

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
    Exp.Channels = [ 1 0 1 0 1  0 1 0 1 0];
    Exp.Inc =  [ 0 0 +1 0 0];
    Exp.t =  [tau1 tau1 T0 tau2 tau2]*1e-3;
    Exp.tp = [16 32 16 16 32]*1e-3;
  end
end

nPoints = 256; dt = 0.08;

nIntervals = numel(Exp.t);

% Compute all CTPs which are detectable
N = nIntervals-1;
[idx{1:N}] = ndgrid(1:4);
CTP = fliplr(reshape(cat(ndims(idx{1}),idx{:}),[],N)) ;
CTP(:,nIntervals) = 4;
nCTP = 4^N;
neCTP = 4^N - 3^N;

co = [0 0 +1 -1];
IncDims = Exp.Inc;
IntervalTimes = Exp.t;
ElCoh = co(CTP);

% Identify all terms that never refocus
% -> echo if total time in - equals total time in +
for c = 1:nCTP

  % Determine number of order changes |delta p|=2
  cc = ElCoh(c,:);
  cc(cc==0) = [];
  Echoes(c) = sum(abs(diff(cc))==2);

  % Determine echo phase contribution for constant intervals
  t0(c) = ElCoh(c,:)*IntervalTimes.';

  % Determine echo phase contribution for sweep intervals
  for d = 1:max(abs(IncDims))
    incdec(c,d) = sum(ElCoh(c,abs(IncDims)==d));
  end

  Refocus(c) = (Echoes(c)>0) && (t0(c)==0) && all(incdec(c,:)==0);

  Crosses(c) = (Echoes(c)>0) && any((t0(c)+nPoints*dt*incdec(c,:))*t0(c)<=0);

end

SelectedCTP = CTP(Refocus,:);


if (Display)
  % Summary
  fprintf('%d pulses\n%d eCTPs\n%d eCTPs give at least one echo\n',...
    nIntervals,nCTP,neCTP);
  fprintf('%d eCTPs lead to an echo always in the detection window\n',sum(Refocus));
  fprintf('%d crossing echoes\n',sum(Crosses)-sum(Refocus));
  % Symbolic representation of CTPs with +, -, a and b
  Str = 'ab+-';
  fprintf('CTP   #echos  time    walk dirs\n');
  for iCTP = 1:nCTP
    if Echoes(iCTP)>0 && Crosses(iCTP)
      fprintf('%s    %d    %+5.5g     ',Str(CTP(iCTP,:)),Echoes(iCTP),t0(iCTP));
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

% Some phase cyclig stuff

%ElCoh = co(CTP);
%deltap = diff(ElCoh,2);
