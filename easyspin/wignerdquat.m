% wignerdquat  Compute Wigner D-matrix element from a given quaternion
%
%   D = wignerdquat(J,m1,m2,q);
%   D = wignerdquat(J,m1,m2,Q);
%
%   Input:
%      j              rank of Wigner matrix
%      q              quaternion, either in 4-component or Pauli matrix
%                     form
%
%   Output:
%      D              Wigner D-matrix

% Implementation based on 
% - Hanson, A.J., Visualizing Quaternions (2006)


function D = wignerdquat(j,m1,m2,Q)
% Calculate the D-matrix using a quaternion polynomial, see Eq. 27.9 in
% reference, where m1=m' and m2=m

if numel(j)~=1 || ~isreal(j) || mod(j,1.0) || (j<0) || ~isnumeric(j)
  error('Index j must be an integer greater than zero.')
end

if numel(m1)~=1 || ~isreal(m1) || mod(m1,1.0) || (abs(m1)>j) || ~isnumeric(m1)
  error('Index m1 must be in the range -j<=m1<=j.')
end

if numel(m2)~=1 || ~isreal(m2) || mod(m2,1.0) || (abs(m2)>j) || ~isnumeric(m2)
  error('Index m2 must be in the range -j<=m2<=j.')
end

if ~isnumeric(Q) || size(Q,1)~=2 || size(Q,2)~=2
  error('Q must be a 2x2x... array.')
end

q = squeeze([Q(1,1,:); Q(1,2,:); Q(2,1,:); Q(2,2,:)]);

% Pre-compute the factorials
prefactor = sqrt(factorial(j+m1)*factorial(j-m1)*factorial(j+m2)*factorial(j-m2));
s = sIndices(j,m1,m2);
powers = repmat(permute([j+m2-s;m1-m2+s;s;j-m1-s],[1,3,2]),[1,size(Q,3),1]);
denoms = repmat(permute([factorial(j+m2-s); factorial(m1-m2+s); factorial(s); factorial(j-m1-s)],[1,3,2]),[1,size(Q,3),1]);

% Manipulate q using copies so that it can be acted upon by powers and denoms
q = repmat(q, [1,1,numel(s)]);

D = prefactor*sum(prod(bsxfun(@rdivide,bsxfun(@power,q,powers),denoms),1),3);

function s = sIndices(j,m1,m2)
% Compute the values for m1, m2, and s for a given j
s = max(0,m2-m1):min(j+m2,j-m1);

end

end
