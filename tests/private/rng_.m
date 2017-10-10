%  rng_(seed, generator) seed the random number generator in a way that is
%                        compatible with old (pre-R2011b) and new Matlab
%                        versions
%
%   seed          non-negative integer
%
%   generator     character array
%
%                 'twister'               Mersenne Twister  (default)
%                 'simdTwister'           SIMD-oriented Fast Mersenne Twister
%                 'combRecursive'         Combined Multiple Recursive
%                 'multFibonacci'         Multiplicative Lagged Fibonacci
%                 'v5uniform'             Legacy MATLAB 5.0 uniform generator
%                 'v5normal'              Legacy MATLAB 5.0 normal generator
%                  'v4'                    Legacy MATLAB 4.0 generator

function rng_(seed, generator)

if seed<0||floor(seed)~=seed
  error('seed must be a non-negative integer.')
end

if exist('generator','var')~=1
  generator = 'twister';
elseif ~ischar(generator)
  error('generator must be a character array.')
end

if verLessThan('matlab','7.12')  % 7.12 = R2011b
  rand(generator,seed)
else
  rng(seed,generator)
end

end

