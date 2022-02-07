function ok = test()

% Test diamond space group (Fd-3m, 227)
R = runprivate('sitetransforms','Fd-3m');

% Combine all rotation matrices as column vectors in a matrix, and sort
for k = 1:24
  RR(k,:) = R{k}(:);
end
RR = sortrows(RR).';

% Reference rotation matrices in same format
RR_ref = [
    -1    -1    -1    -1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     1     1     1
     0     0     0     0    -1    -1    -1    -1     0     0     0     0     0     0     0     0     1     1     1     1     0     0     0     0
     0     0     0     0     0     0     0     0    -1    -1    -1    -1     1     1     1     1     0     0     0     0     0     0     0     0
     0     0     0     0    -1     0     0     1    -1     0     0     1    -1     0     0     1    -1     0     0     1     0     0     0     0
    -1     0     0     1     0     0     0     0     0    -1     1     0     0    -1     1     0     0     0     0     0    -1     0     0     1
     0    -1     1     0     0    -1     1     0     0     0     0     0     0     0     0     0     0    -1     1     0     0    -1     1     0
     0     0     0     0     0     1    -1     0     0    -1     1     0     0     1    -1     0     0    -1     1     0     0     0     0     0
     0    -1     1     0     0     0     0     0     1     0     0    -1    -1     0     0     1     0     0     0     0     0     1    -1     0
     1     0     0    -1    -1     0     0     1     0     0     0     0     0     0     0     0     1     0     0    -1    -1     0     0     1
     ];

ok = all(RR-RR_ref==0);
