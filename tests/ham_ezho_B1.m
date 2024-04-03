function ok = test()

% Assert that ham_ezho with 1st order in B gives the same result as ham_ez.

rng(5);

B = rand(1,3);  % magnetic field, mT

S = 5/2;
Ham110 = rand; % (lB,lS,l)=(1,1,0), m=0; MHz/mT
Ham112 = rand(5,1); % (lB=1,lS,l)=(1,1,2), m= l,...,-l; MHz/mT

% Construct g matrix equivalent to Ham110 and Ham112 input to ham_ezho
g = zeros(3);
g(1,2) = Ham112(5)/sqrt(2);
g(1,3) = Ham112(2)/sqrt(2);
g(2,3) = Ham112(4)/sqrt(2);
g = g + g';
g(3,3) = (-Ham110+sqrt(2)*Ham112(3))/sqrt(3);
g(2,2) = (-sqrt(2)*Ham110-Ham112(3)-sqrt(3)*Ham112(1))/sqrt(6);
g(1,1) = (-sqrt(2)*Ham110-Ham112(3)+sqrt(3)*Ham112(1))/sqrt(6);
g = g *(planck*1e9)/bmagn; % MHz/mT -> unitless

% Assemble spin systems
Sys.S = S;
Sys.Ham110 = Ham110;
Sys.Ham112 = Ham112;

Sys2.S = S;
Sys2.g = g;

% Evaluate using ham_ez and ham_ezho, with and without field input
zHo = ham_ezho(Sys);
H_ezho{1} = -zHo{2}{1};
H_ezho{2} = -zHo{2}{2};
H_ezho{3} = -zHo{2}{3};
H_ezho{4} = ham_ezho(Sys,B);

[H_ez{1},H_ez{2},H_ez{3}] = ham_ez(Sys2);
H_ez{4} = ham_ez(Sys2,B);

% Comparison
threshold = 1e-8;
for n = 4:-1:1
  ok(n) = areequal(H_ezho{n},H_ez{n},threshold,'abs');
end
