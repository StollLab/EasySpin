% stochtraj_wignerdquat  Compute Wigner D-matrix element from a quaternion 
%                        or set of quaternions.
%
%   D = stochtraj_wignerdquat(L,M,K,q);
%
%   Input:
%     L              int
%                    rank of Wigner matrix
%
%     M              int
%
%     K              int
%
%     q              numeric, size = (4,nTraj)
%                    quaternions
%
%
%   Output:
%     Del            complex
%                    Wigner D-matrix element

%   References
%   ----------
%   [1] Hanson, A. J., Ch. 27.2.1, Visualizing Quaternions (2006)


function Del = stochtraj_wignerdquat(L, M, K, q)
% Calculate the D-matrix using a quaternion polynomial, see Eq. 27.9 in
% reference, where M=m' and K=m

if nargin==0, help(mfilename); return; end

persistent cache;


% Error check
% -------------------------------------------------------------------------

if numel(L)~=1 || ~isreal(L) || mod(L,1.0) || (L<0) || ~isnumeric(L)
  error('Index L must be an integer greater than zero.')
end

if numel(M)~=1 || ~isreal(M) || mod(M,1.0) || (abs(M)>L) || ~isnumeric(M)
  error('Index M must be in the range -L<=M<=L.')
end

if numel(K)~=1 || ~isreal(K) || mod(K,1.0) || (abs(K)>L) || ~isnumeric(K)
  error('Index K must be in the range -L<=K<=L.')
end

if ~isnumeric(q) || size(q,1)~=4
  error('q must be an array of size (4,...) array.')
end

% % Adjust to match output of wignerd
% switch ndims(q)
%   case 2, q(1,:) = -q(1,:);
%   case 3, q(1,:,:) = -q(1,:,:);
% end


% Setup
% -------------------------------------------------------------------------

if isempty(cache) || ~all(cache.LMK==[L,M,K]) || ~all(cache.sizeq==size(q))
  cache.sizeq = size(q);
  cache.LMK = [L,M,K];
  cache.prefactor = sqrt(factorial(L+M)*factorial(L-M)*factorial(L+K)*factorial(L-K));
  cache.s = max(0,K-M):min(L+K,L-M);
  cache.powers = repmat(permute([ L+K-cache.s; ...
                                  M-K+cache.s; ...
                                      cache.s; ...
                                  L-M-cache.s ], ...
                          [1,3,2]),[1,size(q,2),1]);
  cache.denoms = repmat(permute([ factorial(L+K-cache.s); ...
                                  factorial(M-K+cache.s); ...
                                      factorial(cache.s); ...
                                  factorial(L-M-cache.s) ], ...
                                 [1,3,2]),[1,size(q,2),1]);
end

% Manipulate q so that it can be acted upon by powers and denoms
Q = [    q(1,:) -  1i*q(4,:); ...
     -1i*q(2,:) -     q(3,:); ...
     -1i*q(2,:) +     q(3,:); ...
         q(1,:) +  1i*q(4,:)];
Q = repmat(Q, [1,1,numel(cache.s)]);

% Calculate
Del = cache.prefactor*sum(...
                        prod(...
                             bsxfun(@rdivide,...
                                    bsxfun(@power,...
                                           Q,cache.powers),...
                                    cache.denoms),...
                             1),...
                        3);
