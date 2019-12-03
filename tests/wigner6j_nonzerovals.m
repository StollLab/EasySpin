function [err,data] = test(opt,olddata)

% various non-zero values
%======================================================
Jj(1,:) = [3 3 2 3/2 1/2 5/2];  val(1) = +sqrt(3/35)/2;
Jj(2,:) = [5/2,2,1/2,2,5/2,1];  val(2) = +sqrt(7)/15;
Jj(3,:) = [5/2,3,1/2,3,5/2,1];  val(3) = -sqrt(10)/21;
Jj(4,:) = [4,7/2,3/2,3/2,1,3];  val(4) = -1/3*sqrt(1/7);
a = 5; b = 3; c = 2;
Jj(5,:) = [a,b,c,0,c,b];        val(5) = (-1)^(a+b+c)/sqrt((2*b+1)*(2*c+1));
Jj(6,:) = [2 2 1 2 1 2];        val(6) = sqrt(7/3)/10;
Jj(7,:) = [6 6 4 9/2 7/2 11/2]; val(7) = 1363/36036;
Jj(8,:) = [8 8 8 8 7 7];        val(8) = sqrt(7/3)*2557/193154;
Jj(9,:) = [1 2 3 4 5 6];        val(9) = sqrt(2/715)/3;

for k = 1:numel(val)
  a(k) = wigner6j(Jj(k,:));
end

%[a;val].'

err = ~areequal(a,val,1e-12,'abs');

data = [];
