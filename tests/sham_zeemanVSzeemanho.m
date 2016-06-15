function [err,data] = test(opt,olddata)
%==================================================================
% Test whether the 1th order in B of zeemanho is identical to the usual
% Zeeman Hamiltonian
%==================================================================

Sys1.S =randi(10)/2;

B = rand(1,3);

Sys1.ZB11.l = [0,2];
Sys1.ZB11.vals{1} = rand; % lB = 1, lS =1 l, = 0, m= 0
Sys1.ZB11.vals{2} = rand(5,1); % lB = 1, lS =1 l, = 2, m= l,...,-l
a = Sys1.ZB11.vals{1};
b = Sys1.ZB11.vals{2};

Sys2.S = Sys1.S;
Sys2.g = zeros(3);

Sys2.g(1,2) = gfree*b(5)/sqrt(2);
Sys2.g(1,3) = gfree*b(2)/sqrt(2);
Sys2.g(2,3) = gfree*b(4)/sqrt(2);

Sys2.g = Sys2.g +Sys2.g'; %symmetrize before diagonal elements

Sys2.g(3,3) = gfree*(-a+sqrt(2)*b(3))/sqrt(3);
Sys2.g(2,2) = gfree*(-sqrt(2)*a-b(3)-sqrt(3)*b(1))/sqrt(6);
Sys2.g(1,1) = gfree*(-sqrt(2)*a-b(3)+sqrt(3)*b(1))/sqrt(6);

Hz{1} = sham(Sys1,B);
H{1} = sham(Sys2,B);

[Hz{2},Hz{3}] = sham(Sys1,B);
[H{2},H{3}] = sham(Sys2,B);

[Hz{4},Hz{5},Hz{6},Hz{7}] = sham(Sys1);
[H{4},H{5},H{6},H{7}] = sham(Sys2);

% test
threshold = 1e-8;
for n = 7:-1:1
  errl(n) = ~areequal(H{n},Hz{n},threshold);
end
err = any(errl);
data =[];


