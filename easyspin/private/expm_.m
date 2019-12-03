% EXPM_    Matrix exponential, using multiple algorithms
%
%  C = expm_(A)
%  C = expm_(A,algo)
%
%  Input:
%    A     matrix to exponentiate
%    algo  algorithm to use:
%          'orig'   Matlab built-in (default)
%          'eig'    using eigenvalues and eigenvectors
%          'fast'   a stripped-down Matlab's built-in, faster than the original
%
%  Output:
%    C     matrix eponential expm(A)

function C = expm_(A,algo)

if nargin<2, algo = 'orig'; end

switch algo
  case 'orig'
    C = expm(A);
  case 'fast'
    C = expm_fast1(A);
  case 'eig'
    [V,D] = eig(A);
    C = V*diag(exp(diag(D)))*V';
end

end

function F = expm_fast1(A)
%EXPM  Matrix exponential.
%   EXPM(A) is the matrix exponential of A and is computed using
%   a scaling and squaring algorithm with a Pade approximation.
%
%   Although it is not computed this way, if A has a full set
%   of eigenvectors V with corresponding eigenvalues D then
%   [V,D] = EIG(A) and EXPM(A) = V*diag(exp(diag(D)))/V.
%
%   EXP(A) computes the exponential of A element-by-element.
%
%   See also LOGM, SQRTM, FUNM.

%   References:
%   N. J. Higham, The scaling and squaring method for the matrix
%      exponential revisited. SIAM J. Matrix Anal. Appl., 26(4), (2005),
%      pp. 1179-1193.
%   A. H. Al-Mohy and N. J. Higham, A new scaling and squaring algorithm
%      for the matrix exponential, SIAM J. Matrix Anal. Appl., 31(3),
%      (2009), pp. 970-989.
%
%   Nicholas J. Higham and Samuel D. Relton
%   Copyright 1984-2015 The MathWorks, Inc.

T = A;

if isdiag(T) % Check if T is diagonal.
  d = diag(T);
  F = diag(exp(full(d)));
  return;
end

% Compute exponential
% Get scaling and Pade parameters.
[s, m, T2, T4, T6] = expm_params(T);

% Rescale the powers of T appropriately.
if s ~= 0
  T = T/(2.^s);
  T2 = T2/2^(s*2);
  T4 = T4/2^(s*4);
  T6 = T6/2^(s*6);
end

% Evaluate the Pade approximant.
switch m
  case 3
    c = [120, 60, 12, 1];
  case 5
    c = [30240, 15120, 3360, 420, 30, 1];
  case 7
    c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
  case 9
    c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
      2162160, 110880, 3960, 90, 1];
  case 13
    c = [64764752532480000, 32382376266240000, 7771770303897600, ...
      1187353796428800,  129060195264000,   10559470521600, ...
      670442572800,      33522128640,       1323241920,...
      40840800,          960960,            16380,  182,  1];
end
I = eye(length(T));
switch m
  case {3, 5, 7, 9}
    Tpowers = {[],T2,[],T4,[],T6};
    strt = length(Tpowers) + 2;
    for k = strt:2:m-1
      Tpowers{k} = Tpowers{k-2}*Tpowers{2};
    end
    U = c(2)*I;
    V = c(1)*I;
    for j = m:-2:3
      U = U + c(j+1)*Tpowers{j-1};
      V = V + c(j)*Tpowers{j-1};
    end
    U = T*U;
  case 13
    U = T * (T6*(c(14)*T6 + c(12)*T4 + c(10)*T2) + c(8)*T6 + c(6)*T4 + c(4)*T2 + c(2)*I);
    V = T6*(c(13)*T6 + c(11)*T4 + c(9)*T2) + c(7)*T6 + c(5)*T4 + c(3)*T2 + c(1)*I;
end
F = (V-U)\(2*U) + I;  %F = (-U+V)\(U+V);

% Squaring phase.
for k = 1:s
  F = F*F;
end
end


function t = ell(T, coeff, m_val)
%ell Function needed to compute optimal parameters.
scaledT = coeff.^(1/(2*m_val+1)) .* abs(T);
alpha = normAm(scaledT,2*m_val+1)/norm(T,1);
t = max(ceil(log2(2*alpha/eps(class(alpha)))/(2*m_val)),0);
end

function [s, m, T2, T4, T6] = expm_params(T)
%expm_params Obtain scaling parameter and order of the Pade approximant.
% Coefficients of backwards error function.
coeff = [1/100800, 1/10059033600, 1/4487938430976000,...
  1/5914384781877411840000, 1/113250775606021113483283660800000000];

s = 0;
% m_val is one of [3 5 7 9 13];
% theta_m for m=1:13.
theta = [%3.650024139523051e-008
  %5.317232856892575e-004
  1.495585217958292e-002  % m_vals = 3
  %8.536352760102745e-002
  2.539398330063230e-001  % m_vals = 5
  %5.414660951208968e-001
  9.504178996162932e-001  % m_vals = 7
  %1.473163964234804e+000
  2.097847961257068e+000  % m_vals = 9
  %2.811644121620263e+000
  %3.602330066265032e+000
  %4.458935413036850e+000
  5.371920351148152e+000];% m_vals = 13

T2 = T*T;
T4 = T2*T2;
T6 = T2*T4;
d4 = norm(T4,1)^(1/4);
d6 = norm(T6,1)^(1/6);
eta1 = max(d4, d6);
if (eta1 <= theta(1) && ell(T, coeff(1), 3) == 0)
  m = 3;
  return;
end
if (eta1 <= theta(2) && ell(T, coeff(2), 5) == 0)
  m = 5;
  return;
end

isSmall = size(T,1) < 150; %Compute matrix power explicitly
if isSmall
  d8 = norm(T4*T4,1)^(1/8);
else
  d8 = normAm(T4, 2)^(1/8);
end
eta3 = max(d6, d8);
if (eta3 <= theta(3) && ell(T, coeff(3), 7) == 0)
  m = 7;
  return;
end
if (eta3 <= theta(4) && ell(T, coeff(4), 9) == 0)
  m = 9;
  return;
end
if isSmall
  d10 = norm(T4*T6,1)^(1/10);
else
  d10 = normAm(T2, 5)^(1/10);
end
eta4 = max(d8, d10);
eta5 = min(eta3, eta4);
s = max(ceil(log2(eta5/theta(5))), 0);
s = s + ell(T/2^s, coeff(5), 13);
if isinf(s)
  % Overflow in ell subfunction. Revert to old estimate.
  [t, s] = log2(norm(T,1)/theta(end));
  s = s - (t == 0.5); % adjust s if normA/theta(end) is a power of 2.
end
m = 13;
end

function [c,mv] = normAm(A,m)
%NORMAM   Estimate of 1-norm of power of matrix.
%   NORMAM(A,M) estimates norm(A^m,1). If A has nonnegative elements
%   the estimate is exact.
%   [C, MV] = NORMAM(A,M) returns the estimate C and the number MV of
%   matrix-vector products computed involving A or A^*.

%   Reference:
%   A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring Algorithm
%      for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
%      970-989, 2009.
%
%   Awad H. Al-Mohy and Nicholas J. Higham
%   Copyright 2014-2015 The MathWorks, Inc.

n = size(A,1);
if n < 50 % Compute matrix power explicitly
  mv = 0;
  c = norm(matlab.internal.math.mpower.viaMtimes(A,m),1);
elseif isreal(A) && all(A(:) >= 0)
  % For positive matrices only.
  e = ones(n,1,class(A));
  for j=1:m
    e = A'*e;
  end
  c = norm(e,inf);
  mv = m;
else
  [c,~,~,it] = normest1(@afun_power);
  mv = it(2)*2*m; % Since t = 2.
end
% End of normAm

  function Z = afun_power(flag,X)
    %afun_power  Function to evaluate matrix products needed by normest1.
    if isequal(flag,'dim')
      Z = n;
    elseif isequal(flag,'real')
      Z = isreal(A);
    else
      if isequal(flag,'notransp')
        for i = 1:m
          X = A*X;
        end
      elseif isequal(flag,'transp')
        for i = 1:m
          X = A'*X;
        end
      end
      Z = X;
    end
  end
end
