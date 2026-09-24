function ok = test()

% Recursion method ('r') agrees with Racah formula ('f') for small j

jmax = 6;
maxerr = 0;
for j1 = 0:1/2:jmax
  for j2 = 0:1/2:jmax
    for j3 = abs(j1-j2):min(j1+j2,jmax)
      for m1 = -j1:j1
        for m2 = -j2:j2
          m3 = -m1-m2;
          if abs(m3)>j3, continue; end
          a = wigner3j(j1,j2,j3,m1,m2,m3,'r');
          b = wigner3j(j1,j2,j3,m1,m2,m3,'f');
          maxerr = max(maxerr,abs(a-b));
        end
      end
    end
  end
end

ok = maxerr<1e-12;

end
