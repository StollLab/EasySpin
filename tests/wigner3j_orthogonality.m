function ok = test()

% Orthogonality relation of 3j symbols (sum over m1 and m2):
%   sum_{m1,m2} (2*j3+1) (j1 j2 j3; m1 m2 m3) (j1 j2 j3'; m1 m2 m3) = delta(j3,j3')
% for all j1, j2 up to 6, integer and half-integer

jmax = 6;
maxerr = 0;
for j1 = 0:1/2:jmax
  for j2 = 0:1/2:jmax
    j3list = abs(j1-j2):j1+j2;
    for m3 = -max(j3list):max(j3list)
      % 3j symbols for all j3 and all (m1,m2) with m1+m2+m3=0
      m1list = -j1:j1;
      T = zeros(numel(j3list),numel(m1list));
      for a = 1:numel(j3list)
        for b = 1:numel(m1list)
          T(a,b) = wigner3j(j1,j2,j3list(a),m1list(b),-m1list(b)-m3,m3);
        end
      end
      % only j3 with |m3|<=j3 contribute
      valid = abs(m3)<=j3list;
      S = diag(2*j3list(valid)+1)*T(valid,:)*T(valid,:).';
      if ~isempty(S)
        maxerr = max(maxerr,max(max(abs(S-eye(nnz(valid))))));
      end
    end
  end
end

ok = maxerr<1e-13;

end
