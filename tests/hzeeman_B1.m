function [err,data] = test(opt,olddata)
% first order in B of highZeeman has to be identical to zeeman

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

Hz = hzeeman(Sys,B);
H = sham(Sys,B);

err = round(sum(sum(Hz-H)),8);
data =[];
end

