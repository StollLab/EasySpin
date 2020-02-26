function ok = test()

tri = sphtri('D2h',2);
ok = all(tri==[1 2 3]);
