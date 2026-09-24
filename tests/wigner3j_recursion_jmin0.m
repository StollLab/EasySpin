function ok = test()

% Recursion method for j2==j3 and m1==0, where the recursion starts at j1=0

j = 30;
m = 7;
a = wigner3j(0,j,j,0,m,-m,'r');
b = (-1)^(j-m)/sqrt(2*j+1);
c = wigner3j(1,j,j,0,m,-m,'r');
d = (-1)^(j-m)*m/sqrt(j*(j+1)*(2*j+1));

% larger j1 (j1 is recursed over, since it is largest)
e = wigner3j(50,j,j,0,m,-m,'r');
f = wigner3j(50,j,j,0,m,-m,'r+'); % same, fast paths don't apply

ok = areequal([a c],[b d],1e-14,'abs') && e==f && e~=0;

end
