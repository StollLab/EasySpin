%  grstat  Calculate the Gelman-Rubin R statistic for a set of Markov
%           chains.
%
%  gr = grstat(A);
%
%  Input:
%     A        NxM array, M Markov chains each with length N
%           or LxNxM array, L instances of M Markov chains with length N
%
%  Output:
%     gr       double, Gelman-Rubin R statistic
%           or Lx1 array, Gelman-Rubin R statistic for each of the L
%              instances

% Implementation based on
%   Gelman and Rubin, Journal of Computational and Graphical Statistics 7,
%   434 (1998)
%     http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.55.1675

function stat = stochtraj_grstat(A)

if (nargin==0), help(mfilename); return; end

ndimsA = ndims(A);
sizeA = size(A);

if ~isnumeric(A)
  error('Must provide a numeric array for input.')
end

if ndimsA>3
  error('Array dimensions greater than 3 are not supported.')
end

if isrow(A) || (ndimsA==3 && ~sizeA(2)>1)
  error('Markov chain length must be greater than one.')
end

if iscolumn(A) || (ndimsA==3 && ~sizeA(3)>1)
  error('Must have more than one Markov chain.')
end

if ndimsA==2
  % Scalar Markov chains
  stat = find_gr(A);
elseif ndimsA==3
  % Vector Markov chains
  nComps = size(A,1);
  stat = zeros(nComps,1);
  for iComp=1:nComps
    stat(iComp) = find_gr(squeeze(A(iComp,:,:)));
  end
else
  error('Number of input array dimensions must be 2 or 3.')
end

  function gr = find_gr(X)
    
    n = size(X,1);  % chain length
    m = size(X,2);  % number of chains
    
    % Inter-chain variance
    B = n*var(mean(X,1));  % B from reference
    
    % Intra-chain variance
    W = mean(var(X,0,1));  % W from reference
    
    sigmasq = (n-1)/n*W + B/n;
    Vhat = sigmasq + B/m/n;
    % Vhat = (n-1)/n*W + (m+1)/(m*n)*B;
    
    % Calculate the R statistic
    gr = Vhat/W;  %
    
  end

end