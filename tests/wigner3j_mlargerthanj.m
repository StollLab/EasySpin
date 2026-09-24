function ok = test()

% 3j symbols with |m|>j are zero

a = [wigner3j(30,30,30,31,-15,-16) ...
     wigner3j(30,30,30,-15,31,-16) ...
     wigner3j(3,3,3,4,-2,-2)];

ok = all(a==0);

end
