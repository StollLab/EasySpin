function ok = test()

% Recursion method for half-integer j, compared with Racah formula ('f'),
% which is accurate for j up to about 20

a = [wigner3j(25/2,21/2,7,-3/2,1/2,1,'r') ...
     wigner3j(35/2,19/2,13,5/2,-9/2,2,'r') ...
     wigner3j(13,39/2,31/2,-4,15/2,-7/2,'r')];
b = [wigner3j(25/2,21/2,7,-3/2,1/2,1,'f') ...
     wigner3j(35/2,19/2,13,5/2,-9/2,2,'f') ...
     wigner3j(13,39/2,31/2,-4,15/2,-7/2,'f')];

ok = areequal(a,b,1e-12,'abs') && all(a~=0);

end
