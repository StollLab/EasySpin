function ok = test()

% Test for some randomly selected parameter combinations

a(1) = wigner3j(11/2,5,9/2,7/2,-3,-1/2);
b(1) = 3/2/sqrt(286);
a(2) = wigner3j(2,2,2,1,0,-1);
b(2) = 1/sqrt(70);
a(3) = wigner3j(7,5,4,0,0,0);
b(3) = 2/3*sqrt(70/2431);
a(4) = wigner3j(5,5,5,2,1,-3);
b(4) = -sqrt(5/429);
a(5) = wigner3j(10,10,10,0,0,0);
b(5) = -126*sqrt(7/33393355);
a(6) = wigner3j(20,20,20,0,0,0);
b(6) = 92378*sqrt(13/126902256462843);
a(7) = wigner3j(10,10,9,1,0,-1);
b(7) = 21*sqrt(7/1077205);
a(8) = wigner3j(4,4,4,4,0,-4);
b(8) = sqrt(14/143)/3;
a(9) = wigner3j(4,4,4,3,1,-4);
b(9) = -sqrt(35/143)/3;

ok = areequal(a,b,1e-6,'abs');
