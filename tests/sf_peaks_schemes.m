function ok = test() %#ok<FNDEF>

% sf_peaks bins the peaks of each incrementation scheme like the direct sum
% over all state indices: the amplitude is the product of the matrix
% elements along the pathway (see the sf_peaks header), and the bin is
% fmod(-nu*dt,1)*nPoints, rounded toward zero, with negative bins wrapped
% around.
% - random complex, non-Hermitian matrices, so that every matrix element
%   must enter with the right indices (no symmetry hides an error)
% - random manifolds (alpha/beta) per incremented delay, buffer sizes that
%   are not powers of 2
% - non-finite energies or time increments and empty buffers are rejected

rng(1);
n = 3;  % number of nuclear states
dt = [0.008 0.013];  % µs
nPoints = [210 150];
schemes = [1 2 3 11 12 14 15 17];
nFree = [1 2 2 2 3 3 4 4];  % number of incremented delays

ok = false(1,numel(schemes));
for s = 1:numel(schemes)
  scheme = schemes(s);
  nf = nFree(s);
  Ea = 10*randn(n,1);  % MHz
  Eb = 10*randn(n,1);  % MHz
  iL = randi(2,1,nf);
  iR = randi(2,1,nf);
  G = randn(n)+1i*randn(n);
  D = randn(n)+1i*randn(n);
  ML = cell(1,nf-1);
  MR = cell(1,nf-1);
  for k = 1:nf-1
    ML{k} = randn(n)+1i*randn(n);
    MR{k} = randn(n)+1i*randn(n);
  end

  if scheme>10
    d = dt;
    siz = nPoints;
  else
    d = dt(1);
    siz = [1 nPoints(1)];
  end
  bufferRe = zeros(siz);
  bufferIm = zeros(siz);
  runprivate('sf_peaks',scheme,bufferRe,bufferIm,d,iL,iR,Ea,Eb,G,D,ML{:},MR{:});
  B = complex(bufferRe,bufferIm);

  [amp,nu1,nu2] = peaks(scheme,Ea,Eb,iL,iR,G,D,ML,MR);
  b1 = bins(nu1,dt(1),siz(1+(scheme<10)));
  if scheme>10
    b2 = bins(nu2,dt(2),siz(2));
    Bref = accumarray([b1 b2]+1,amp,siz);
  else
    Bref = accumarray(b1+1,amp,[siz(2) 1]).';
  end
  ok(s) = areequal(B,Bref,1e-12,'rel');
end

% very large phases (|nu*dt| >= 2^52 cycles, integers): all peaks in bin 0
EaL = 1e20*[1; 2; 3];  % MHz
EbL = 1e20*[3; 1; 2];  % MHz
bufferRe = zeros(1,10);
bufferIm = zeros(1,10);
runprivate('sf_peaks',1,bufferRe,bufferIm,dt(1),1,2,EaL,EbL,G,D);
[amp,nu1] = peaks(1,EaL,EbL,1,2,G,D,{},{});
Bref = accumarray(bins(nu1,dt(1),10)+1,amp,[10 1]).';
ok(end+1) = areequal(complex(bufferRe,bufferIm),Bref,1e-12,'rel');

% invalid inputs: non-finite energies or time increment, empty buffer
Ea = [1; 2; 3];  % MHz
Eb = [-1; 0; 4];  % MHz
args = {{[1; NaN; 3],Eb,dt(1),zeros(1,10)}, {Ea,[-1; Inf; 4],dt(1),zeros(1,10)}, ...
  {Ea,Eb,NaN,zeros(1,10)}, {Ea,Eb,dt(1),zeros(1,0)}};
messages = {'must be finite','must be finite','must be finite','must not be empty'};
for k = 1:numel(args)
  [Ea_,Eb_,dt_,buffer] = args{k}{:};
  try
    runprivate('sf_peaks',1,buffer,zeros(size(buffer)),dt_,1,1,Ea_,Eb_,G,D);
    ok(end+1) = false; %#ok<AGROW>
  catch err
    ok(end+1) = contains(err.message,messages{k}); %#ok<AGROW>
  end
end

end

%-------------------------------------------------------------------------------
% bins: fmod(-nu*dt,1)*nPoints, rounded toward zero, negative bins wrapped
function b = bins(nu,dt,nPoints)
b = mod(fix(rem(-nu*dt,1)*nPoints),nPoints);
end

%-------------------------------------------------------------------------------
% peaks: amplitudes and frequencies (dimensions 1 and 2) of all peaks, as
% column vectors over all state indices
function [amp,nu1,nu2] = peaks(scheme,Ea,Eb,iL,iR,G,D,ML,MR)
n = numel(Ea);
E = {Ea,Eb};
El = @(k) E{iL(k)};  % energies of the left and right manifolds of delay k
Er = @(k) E{iR(k)};
el = @(A,r,c) A(sub2ind([n n],r,c));  % elements A(r,c) for index arrays
nu2 = [];
switch scheme
  case 1
    [j,k] = ndgrid(1:n);
    amp = el(D,k,j).*el(G,j,k);
    nu1 = col(El(1),j)-col(Er(1),k);
  case {2,3,11}
    [i,j,k,l] = ndgrid(1:n);
    amp = el(ML{1},i,j).*el(G,j,k).*el(MR{1},k,l).*el(D,l,i);
    nuA = col(El(1),j)-col(Er(1),k);  % first delay
    nuB = col(El(2),i)-col(Er(2),l);  % second delay
    switch scheme
      case 2, nu1 = nuB + nuA;
      case 3, nu1 = -col(El(2),i)+col(El(1),j)-col(Er(1),k)+col(Er(2),l);  % nuA-nuB
      case 11, nu1 = nuA; nu2 = nuB;
    end
  case {12,14}
    [i,j,k,l,m,nn] = ndgrid(1:n);
    amp = el(D,nn,i).*el(ML{2},i,j).*el(ML{1},j,k).*el(G,k,l).*el(MR{1},l,m).*el(MR{2},m,nn);
    nuA = col(El(1),k)-col(Er(1),l);
    nuB = col(El(2),j)-col(Er(2),m);
    nuC = col(El(3),i)-col(Er(3),nn);
    if scheme==12
      nu1 = nuA + nuC; nu2 = nuB;
    else
      nu1 = nuA + nuB; nu2 = nuC;
    end
  case {15,17}
    [i,j,k,l,m,nn,o,p] = ndgrid(1:n);
    amp = el(D,p,i).*el(ML{3},i,j).*el(ML{2},j,k).*el(ML{1},k,l).*el(G,l,m).*...
      el(MR{1},m,nn).*el(MR{2},nn,o).*el(MR{3},o,p);
    nuA = col(El(1),l)-col(Er(1),m);
    nuB = col(El(2),k)-col(Er(2),nn);
    nuC = col(El(3),j)-col(Er(3),o);
    nuD = col(El(4),i)-col(Er(4),p);
    if scheme==15
      nu1 = nuA + nuD; nu2 = nuB + nuC;
    else
      nu1 = nuA + nuB; nu2 = nuC + nuD;
    end
end
amp = amp(:);
nu1 = nu1(:);
nu2 = nu2(:);
end

% elements of a vector for an index array
function v = col(x,idx)
v = x(idx);
end
