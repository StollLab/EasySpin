function ok = test() %#ok<FNDEF>

% sf_peaks bins the peaks of a pathway like the direct sum over all state
% indices: the amplitude is the product of the matrix elements along the
% pathway (see the sf_peaks header), and the frequency along dimension d is
% nu_d = sum_k IncScheme(k,d)*nu_k, with the incrementation scheme
% IncScheme. A peak is binned at fmod(-nu_d*dt(d),1)*nPoints, rounded to
% the nearest integer with ties away from zero, with negative bins wrapped
% around.
% - incrementation schemes of 1 to 4 intervals in 1 or 2 dimensions: the
%   classic ones (3p and 2p ESEEM, HYSCORE, 2D 3p ESEEM, 2D-CP, ...), and
%   schemes with unequal, negative and non-integer steps, intervals
%   incremented in both dimensions (also all of them, so that dimension 2
%   is binned per peak next to tabulated dimension-1 bins), dimension 2
%   first
% - random complex, non-Hermitian matrices, so that every matrix element
%   must enter with the right indices (no symmetry hides an error)
% - diagonal G and transfer matrices (zero amplitudes, skipped early)
% - random manifolds (alpha/beta) per incremented interval, buffer sizes
%   that are not powers of 2
% - exact ties (frequencies on a grid of half bins, exact in binary)
% - real (not complex) matrices
% - non-finite energies, time increments or incrementation schemes, empty
%   or real buffers, manifolds other than 1 and 2, and more than 4
%   intervals are rejected

rng(1);
n = 3;  % number of nuclear states
dt = [0.008 0.013];  % µs
nPoints = [210 150];
schemes = {1, [1;1], [1;-1], [1 0;0 1], [1 0;0 1;1 0], [1 0;1 0;0 1], ...
  [1 0;0 1;0 1;1 0], [1 0;1 0;0 1;0 1], [1;2;1], [1;2;2;1], [1;2.5], ...
  [1 0;0 1;1 1], [0 1;1 0], [1 0;0 1;1 0;0 1], [1 -1;1 1], [1 1], ...
  [0 1;1 1], [0 1;1 1;1 1]};

ok = false(1,numel(schemes));
for s = 1:numel(schemes)
  IncScheme = schemes{s};
  [K,nDims] = size(IncScheme);
  [Ea,Eb,iL,iR,G,D,ML,MR] = randomPathway(n,K);
  if nDims==2
    siz = nPoints;
  else
    siz = [1 nPoints(1)];
  end
  B = complex(zeros(siz));
  runprivate('sf_peaks',B,dt(1:nDims),IncScheme,iL,iR,Ea,Eb,G,D,ML{:},MR{:});
  Bref = directSum(siz,dt(1:nDims),IncScheme,iL,iR,Ea,Eb,G,D,ML,MR);
  ok(s) = areequal(B,Bref,1e-12,'rel');
end

% diagonal G and transfer matrices, as from pulses that keep the electron
% in its manifold: many zero amplitudes, skipped early
for IncScheme = {[1;1], [1 0;0 1;0 1;1 0]}
  [K,nDims] = size(IncScheme{1});
  [Ea,Eb,iL,iR,G,D,ML,MR] = randomPathway(n,K);
  G = diag(diag(G));
  ML = cellfun(@(M)diag(diag(M)),ML,'UniformOutput',false);
  MR = cellfun(@(M)diag(diag(M)),MR,'UniformOutput',false);
  if nDims==2, siz = nPoints; else, siz = [1 nPoints(1)]; end
  B = complex(zeros(siz));
  runprivate('sf_peaks',B,dt(1:nDims),IncScheme{1},iL,iR,Ea,Eb,G,D,ML{:},MR{:});
  Bref = directSum(siz,dt(1:nDims),IncScheme{1},iL,iR,Ea,Eb,G,D,ML,MR);
  ok(end+1) = areequal(B,Bref,1e-12,'rel'); %#ok<AGROW>
end

% very large phases (|nu*dt| >= 2^52 cycles, integers): all peaks in bin 0
[~,~,~,~,G,D] = randomPathway(n,1);
EaL = 1e20*[1; 2; 3];  % MHz
EbL = 1e20*[3; 1; 2];  % MHz
B = complex(zeros(1,10));
runprivate('sf_peaks',B,dt(1),1,1,2,EaL,EbL,G,D);
Bref = directSum([1 10],dt(1),1,1,2,EaL,EbL,G,D,{},{});
ok(end+1) = areequal(B,Bref,1e-12,'rel');

% exact ties: with dt = 1/4 and nPoints = 8, energies that are multiples of
% 1/4 MHz put every peak at an integer or half-integer bin position
[~,~,~,~,G,D,ML,MR] = randomPathway(n,2);
EaT = [0; 0.25; 1.75];  % MHz
EbT = [0.5; -1.25; 3];  % MHz
B = complex(zeros(1,8));
runprivate('sf_peaks',B,0.25,[1;1],[1 2],[2 1],EaT,EbT,G,D,ML{1},MR{1});
Bref = directSum([1 8],0.25,[1;1],[1 2],[2 1],EaT,EbT,G,D,ML,MR);
ok(end+1) = areequal(B,Bref,1e-12,'rel');

% real (not complex) matrices, as with mixed real and complex inputs
[Ea,Eb,~,~,G,D,ML,MR] = randomPathway(n,2);
GR = real(G);
MR1 = real(MR{1});
B = complex(zeros(1,nPoints(1)));
runprivate('sf_peaks',B,dt(1),[1;1],[1 2],[2 1],Ea,Eb,GR,D,ML{1},MR1);
Bref = directSum([1 nPoints(1)],dt(1),[1;1],[1 2],[2 1],Ea,Eb,GR,D,ML,{MR1});
ok(end+1) = areequal(B,Bref,1e-12,'rel');

% invalid inputs: non-finite energies, time increment or incrementation scheme,
% empty or real buffer
Ea = [1; 2; 3];  % MHz
Eb = [-1; 0; 4];  % MHz
buffer = complex(zeros(1,10));
args = {{[1; NaN; 3],Eb,dt(1),1,buffer}, {Ea,[-1; Inf; 4],dt(1),1,buffer}, ...
  {Ea,Eb,NaN,1,buffer}, {Ea,Eb,dt(1),NaN,buffer}, {Ea,Eb,dt(1),1,complex(zeros(1,0))}, ...
  {Ea,Eb,dt(1),1,zeros(1,10)}};
messages = {'must be finite','must be finite','must be finite','must be finite',...
  'must not be empty','must be a complex double array'};
for k = 1:numel(args)
  [Ea_,Eb_,dt_,scheme_,buffer] = args{k}{:};
  try
    runprivate('sf_peaks',buffer,dt_,scheme_,1,1,Ea_,Eb_,G,D);
    ok(end+1) = false; %#ok<AGROW>
  catch err
    ok(end+1) = contains(err.message,messages{k}); %#ok<AGROW>
  end
end

% manifolds other than 1 (alpha) and 2 (beta) are rejected
for idxFree = {{3,1}, {1,0}}
  try
    runprivate('sf_peaks',complex(zeros(1,10)),dt(1),1,idxFree{1}{:},Ea,Eb,G,D);
    ok(end+1) = false; %#ok<AGROW>
  catch err
    ok(end+1) = contains(err.message,'1 (alpha) or 2 (beta)'); %#ok<AGROW>
  end
end

% more than 4 incremented intervals are rejected
try
  runprivate('sf_peaks',complex(zeros(1,10)),dt(1),ones(5,1),ones(1,5),ones(1,5),Ea,Eb,G,D,G,G,G,G,G,G,G,G);
  ok(end+1) = false;
catch err
  ok(end+1) = contains(err.message,'1 to 4 rows');
end

end

%-------------------------------------------------------------------------------
% randomPathway: random energies, manifolds and complex non-Hermitian
% matrices for a pathway with K incremented intervals
function [Ea,Eb,iL,iR,G,D,ML,MR] = randomPathway(n,K)
Ea = 10*randn(n,1);  % MHz
Eb = 10*randn(n,1);  % MHz
iL = randi(2,1,K);
iR = randi(2,1,K);
G = randn(n)+1i*randn(n);
D = randn(n)+1i*randn(n);
ML = cell(1,K-1);
MR = cell(1,K-1);
for k = 1:K-1
  ML{k} = randn(n)+1i*randn(n);
  MR{k} = randn(n)+1i*randn(n);
end
end

%-------------------------------------------------------------------------------
% directSum: buffer from the direct sum over all states a1,b1,...,aK,bK.
% The frequencies of the intervals are added in the loop order of sf_peaks
% (intervals that change dimension 2 first), which matters only for exact
% ties.
function B = directSum(siz,dt,IncScheme,iL,iR,Ea,Eb,G,D,ML,MR)
n = numel(Ea);
[K,nDims] = size(IncScheme);
E = {Ea,Eb};
st = cell(1,2*K);
[st{:}] = ndgrid(1:n);  % st{2k-1} = left state a_k, st{2k} = right state b_k
el = @(A,r,c) A(sub2ind([n n],r,c));  % elements A(r,c) for index arrays
amp = el(G,st{1},st{2}).*el(D,st{2*K},st{2*K-1});
for k = 1:K-1
  amp = amp.*el(ML{k},st{2*k+1},st{2*k-1}).*el(MR{k},st{2*k},st{2*k+2});
end
if nDims==2
  loopOrder = [find(IncScheme(:,2)~=0); find(IncScheme(:,2)==0)].';
else
  loopOrder = 1:K;
end
nu = zeros(numel(amp),nDims);
for k = loopOrder
  nuk = E{iL(k)}(st{2*k-1}(:)) - E{iR(k)}(st{2*k}(:));
  nu = nu + nuk.*IncScheme(k,:);
end
if nDims==1, nBins = prod(siz); else, nBins = siz; end
b = mod(round(rem(-nu.*dt,1).*nBins),nBins);
if nDims==1
  B = accumarray(b+1,amp(:),[nBins 1]).';
else
  B = accumarray(b+1,amp(:),siz);
end
end
