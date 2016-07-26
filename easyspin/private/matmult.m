%  matmult  Perform matrix multiplication across stacks of 3x3 matrices.
%
%  C = matmult(A, B);
%
%  Input:
%     A, B           3x3xnxm arrays
%
%  Output:
%     C              3x3xnxm array

function C = matmult(A, B)

Ashape = size(A);
Bshape = size(B);

if not(isequal(Ashape, Bshape))
  error('Input array shapes do not match.')
end

if not(Ashape(1) == Ashape(2))
  error('The size of the first two dimensions must be equal.')
end

if not(or(Ashape(1) == 2, Ashape(1) == 3))
  error('The size of the first two dimensions must be 2x2 or 3x3.')
end

if or(ndims(A) ~= 4, ndims(B) ~= 4)
    error('Input arrays each must have 4 dimensions.')
end

% Initialize the matrix multiplication output
C = A;

if size(A, 1) == 2
  
  C(1,1,:,:) = A(1,1,:,:).*B(1,1,:,:) + A(1,2,:,:).*B(2,1,:,:);
  C(1,2,:,:) = A(1,1,:,:).*B(1,2,:,:) + A(1,2,:,:).*B(2,2,:,:);

  C(2,1,:,:) = A(2,1,:,:).*B(1,1,:,:) + A(2,2,:,:).*B(2,1,:,:);
  C(2,2,:,:) = A(2,1,:,:).*B(1,2,:,:) + A(2,2,:,:).*B(2,2,:,:);

elseif size(A, 1) == 3

  C(1,1,:,:) = A(1,1,:,:).*B(1,1,:,:) + A(1,2,:,:).*B(2,1,:,:) + A(1,3,:,:).*B(3,1,:,:);
  C(1,2,:,:) = A(1,1,:,:).*B(1,2,:,:) + A(1,2,:,:).*B(2,2,:,:) + A(1,3,:,:).*B(3,2,:,:);
  C(1,3,:,:) = A(1,1,:,:).*B(1,3,:,:) + A(1,2,:,:).*B(2,3,:,:) + A(1,3,:,:).*B(3,3,:,:);

  C(2,1,:,:) = A(2,1,:,:).*B(1,1,:,:) + A(2,2,:,:).*B(2,1,:,:) + A(2,3,:,:).*B(3,1,:,:);
  C(2,2,:,:) = A(2,1,:,:).*B(1,2,:,:) + A(2,2,:,:).*B(2,2,:,:) + A(2,3,:,:).*B(3,2,:,:);
  C(2,3,:,:) = A(2,1,:,:).*B(1,3,:,:) + A(2,2,:,:).*B(2,3,:,:) + A(2,3,:,:).*B(3,3,:,:);

  C(3,1,:,:) = A(3,1,:,:).*B(1,1,:,:) + A(3,2,:,:).*B(2,1,:,:) + A(3,3,:,:).*B(3,1,:,:);
  C(3,2,:,:) = A(3,1,:,:).*B(1,2,:,:) + A(3,2,:,:).*B(2,2,:,:) + A(3,3,:,:).*B(3,2,:,:);
  C(3,3,:,:) = A(3,1,:,:).*B(1,3,:,:) + A(3,2,:,:).*B(2,3,:,:) + A(3,3,:,:).*B(3,3,:,:);

end
% Initialize the matrix multiplication output
% if ndims(A) == 3
%     C = A;
% elseif ndims(B) == 3
%     C = B;
% else
%     error('Neither of the input arrays has 3 dimensions.')
% end
%     
% C(1,1,:) = bsxfun(@times, A(1,1,:), B(1,1,:)) ...
%          + bsxfun(@times, A(1,2,:), B(2,1,:)) ...
%          + bsxfun(@times, A(1,3,:), B(3,1,:));
% C(1,2,:) = bsxfun(@times, A(1,1,:), B(1,2,:)) ...
%          + bsxfun(@times, A(1,2,:), B(2,2,:)) ...
%          + bsxfun(@times, A(1,3,:), B(3,2,:));
% C(1,3,:) = bsxfun(@times, A(1,1,:), B(1,3,:)) ...
%          + bsxfun(@times, A(1,2,:), B(2,3,:)) ...
%          + bsxfun(@times, A(1,3,:), B(3,3,:));
% 
% C(2,1,:) = bsxfun(@times, A(2,1,:), B(1,1,:)) ...
%          + bsxfun(@times, A(2,2,:), B(2,1,:)) ...
%          + bsxfun(@times, A(2,3,:), B(3,1,:));
% C(2,2,:) = bsxfun(@times, A(2,1,:), B(1,2,:)) ...
%          + bsxfun(@times, A(2,2,:), B(2,2,:)) ...
%          + bsxfun(@times, A(2,3,:), B(3,2,:));
% C(2,3,:) = bsxfun(@times, A(2,1,:), B(1,3,:)) ...
%          + bsxfun(@times, A(2,2,:), B(2,3,:)) ...
%          + bsxfun(@times, A(2,3,:), B(3,3,:));
% 
% C(3,1,:) = bsxfun(@times, A(3,1,:), B(1,1,:))
%          + bsxfun(@times, A(3,2,:), B(2,1,:)) ...
%          + bsxfun(@times, A(3,3,:), B(3,1,:));
% C(3,2,:) = bsxfun(@times, A(3,1,:), B(1,2,:)) ...
%          + bsxfun(@times, A(3,2,:), B(2,2,:)) ...
%          + bsxfun(@times, A(3,3,:), B(3,2,:));
% C(3,3,:) = bsxfun(@times, A(3,1,:), B(1,3,:)) ...
%          + bsxfun(@times, A(3,2,:), B(2,3,:)) ...
%          + bsxfun(@times, A(3,3,:), B(3,3,:));

end