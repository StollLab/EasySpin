function ok = test()

rng(36345);

threshold = 1e-10;

ok = true;

a = [0 pi 0];
b = eulang(erot(a));
ok = ok & areequal(a,b,threshold,'abs');

a = [0 0 0];
b = eulang(erot(a));
ok = ok & areequal(a,b,threshold,'abs');

for k = 1:10
  a = rand(1,3)*pi.*[2 1 2];
  b = eulang(erot(a));
  ok = ok && areequal(a,b,threshold,'abs');
end
