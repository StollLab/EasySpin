function [err,data] = test(opt,olddata)
%==================================================================
% Test whether the 1th order in B of zeemanho is identical to the usual
% Zeeman Hamiltonian
%==================================================================

Sys.S =randi(20)/2;

B = rand(1,3);

Sys.ZB11.l = [0,2];
Sys.ZB11.vals{1} = rand; % lB = 1, lS =1 l, = 0, m= 0
Sys.ZB11.vals{2} = rand(5,1); % lB = 1, lS =1 l, = 2, m= l,...,-l
a = Sys.ZB11.vals{1};
b = Sys.ZB11.vals{2};

Sys.g = zeros(3);

Sys.g(1,2) = gfree*b(5)/sqrt(2);
Sys.g(1,3) = gfree*b(2)/sqrt(2);
Sys.g(2,3) = gfree*b(4)/sqrt(2);

Sys.g = Sys.g +Sys.g'; %symmetrize before diagonal elements

Sys.g(3,3) = gfree*(-a+sqrt(2)*b(3))/sqrt(3);
Sys.g(2,2) = gfree*(-sqrt(2)*a-b(3)-sqrt(3)*b(1))/sqrt(6);
Sys.g(1,1) = gfree*(-sqrt(2)*a-b(3)+sqrt(3)*b(1))/sqrt(6);

Hz{1} = zeemanho(Sys,B);
H{1} = zeeman(Sys,B);
zHo = zeemanho(Sys);
for n = 3:-1:1
  Hz{n+1}=zHo{2}{n};
end
[H{2},H{3},H{4}] = zeeman(Sys);

% test
threshold = 1e-8;
for n = 4:-1:1
  errl(n) = ~areequal(H{n},Hz{n},threshold);
end
err = any(errl);
data =[];


