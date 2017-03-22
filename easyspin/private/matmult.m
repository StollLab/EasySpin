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
% Perform 2x2 or 3x3 matrix multiplication across multidimensional arrays

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
  %error('Input arrays each must have 4 dimensions.')
end

% Initialize the matrix multiplication output
C = A;


idx = cell(1, ndims(A));
idx(:) = {':'};

A11 = A(1,1,idx{3:end});
A12 = A(1,2,idx{3:end});
A21 = A(2,1,idx{3:end});
A22 = A(2,2,idx{3:end});

B11 = B(1,1,idx{3:end});
B12 = B(1,2,idx{3:end});
B21 = B(2,1,idx{3:end});
B22 = B(2,2,idx{3:end});

if size(A, 1) == 2
  
  C(1,1,idx{3:end}) = A11.*B11 + A12.*B21;
  C(1,2,idx{3:end}) = A11.*B12 + A12.*B22;

  C(2,1,idx{3:end}) = A21.*B11 + A22.*B21;
  C(2,2,idx{3:end}) = A21.*B12 + A22.*B22;

elseif size(A, 1) == 3

  A13 = A(1,3,idx{3:end});
  A23 = A(2,3,idx{3:end});
  A31 = A(3,1,idx{3:end});
  A32 = A(3,2,idx{3:end});
  A33 = A(3,3,idx{3:end});

  B13 = B(1,3,idx{3:end});
  B23 = B(2,3,idx{3:end});
  B31 = B(3,1,idx{3:end});
  B32 = B(3,2,idx{3:end});
  B33 = B(3,3,idx{3:end});

  C(1,1,idx{3:end}) = A11.*B11 + A12.*B21 + A13.*B31;
  C(1,2,idx{3:end}) = A11.*B12 + A12.*B22 + A13.*B32;
  C(1,3,idx{3:end}) = A11.*B13 + A12.*B23 + A13.*B33;

  C(2,1,idx{3:end}) = A21.*B11 + A22.*B21 + A23.*B31;
  C(2,2,idx{3:end}) = A21.*B12 + A22.*B22 + A23.*B32;
  C(2,3,idx{3:end}) = A21.*B13 + A22.*B23 + A23.*B33;

  C(3,1,idx{3:end}) = A31.*B11 + A32.*B21 + A33.*B31;
  C(3,2,idx{3:end}) = A31.*B12 + A32.*B22 + A33.*B32;
  C(3,3,idx{3:end}) = A31.*B13 + A32.*B23 + A33.*B33;

end



end
