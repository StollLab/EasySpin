function ok = test()

% Test that all input ways for tdm orientation are equivalent

ori = [10 30 120]*pi/180;
k = [78 211]*pi/180;
alpha = 17*pi/180;

% Three equivalent ways to input the orientation of the tdm
tdm_lettercode = 'xyz';
tdm_n = [1;1;1]/sqrt(3);
[tdm_phi,tdm_theta] = vec2ang(tdm_n);

w_lettercode = photoselect(tdm_lettercode,ori,k,alpha);
w_n = photoselect(tdm_n,ori,k,alpha);
w_angles = photoselect([tdm_phi,tdm_theta],ori,k,alpha);

threshold = 1e-10;
ok(1) = areequal(w_n,w_lettercode,threshold,'rel');
ok(2) = areequal(w_angles,w_lettercode,threshold,'rel');
