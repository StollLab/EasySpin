function ok = test()

% Recursion method for half-integer j, compared with exact values

a = [wigner3j(25/2,21/2,7,-3/2,1/2,1,'') ...
     wigner3j(35/2,19/2,13,5/2,-9/2,2,'') ...
     wigner3j(13,39/2,31/2,-4,15/2,-7/2,'')];
b = [195*sqrt(93501394)/46750697 ...
     9199*sqrt(42815967285)/36699400530 ...
     -11353*sqrt(14640232918773605)/37443050943155];

ok = areequal(a,b,1e-14,'abs');

end
