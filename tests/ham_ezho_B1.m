function ok = test()

% Test whether the 1th order in B of ham_ezho is identical to the usual
% Zeeman Hamiltonian

rng(5);

Sys.S = ceil(rand*20)/2;

B = rand(1,3);

Sys.Ham110 = rand; % lB=1, lS=1, l=0, m= 0
Sys.Ham112 = rand(5,1); % lB=1, lS=1, l=2, m= l,...,-l
a = Sys.Ham110;
b = Sys.Ham112;

Sys2.S = Sys.S;
Sys2.g = zeros(3);

Sys2.g(1,2) = b(5)/sqrt(2);
Sys2.g(1,3) = b(2)/sqrt(2);
Sys2.g(2,3) = b(4)/sqrt(2);

Sys2.g = Sys2.g +Sys2.g'; %symmetrize before diagonal elements

Sys2.g(3,3) = (-a+sqrt(2)*b(3))/sqrt(3);
Sys2.g(2,2) = (-sqrt(2)*a-b(3)-sqrt(3)*b(1))/sqrt(6);
Sys2.g(1,1) = (-sqrt(2)*a-b(3)+sqrt(3)*b(1))/sqrt(6);

% Hamm11l are given in MHz/mT, conversion require division by the Bohr magneton
Sys2.g = Sys2.g *(planck*1e9)/bmagn;

Hz{1} = ham_ezho(Sys,B);
zHo = ham_ezho(Sys);
Hz{2} = -zHo{2}{1};
Hz{3} = -zHo{2}{2};
Hz{4} = -zHo{2}{3};

H{1} = ham_ez(Sys2,B);
[H{2},H{3},H{4}] = ham_ez(Sys2);

% test
threshold = 1e-8;
for n = 4:-1:1
  ok(n) = areequal(H{n},Hz{n},threshold,'abs');
end
