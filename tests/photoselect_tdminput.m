function ok = test()

% Test that all input ways for tdm orientation are equivalent

ori = [10 30 120]*pi/180;
k = [78 211]*pi/180;
alpha = 17*pi/180;

lettercode = 'xyz';
n = [1;1;1]/sqrt(3);
[phi,theta] = vec2ang(n);

w1 = photoselect(lettercode,ori,k,alpha);
w2 = photoselect(n,ori,k,alpha);
w3 = photoselect([phi,theta],ori,k,alpha);

threshold = 1e-10;
ok(1) = areequal(w2,w1,threshold,'rel');
ok(2) = areequal(w3,w1,threshold,'rel');
