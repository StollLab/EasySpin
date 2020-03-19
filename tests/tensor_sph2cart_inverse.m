function ok = test()

T = [1 2 3; 4 5 6; 7 8 9];

[T0,T1,T2] = tensor_cart2sph(T);

T_ = tensor_sph2cart(T0,T1,T2);

ok = areequal(T_,T,1e-14,'rel');
