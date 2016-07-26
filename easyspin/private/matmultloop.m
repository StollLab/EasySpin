function C = matmultloop(A, B)

if not(isequal(size(A), size(B)))
  error('Input matrix shapes do not match.') 
end

C = A;

for k = 1:size(A, 4)
  for j = 1:size(A, 3)
    C(:, :, j, k) = A(:, :, j, k)*B(:, :, j, k);
  end
end

end