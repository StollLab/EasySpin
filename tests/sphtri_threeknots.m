function ok = test()

tri = sphtri('D2h',3);
tri0 = [           
           1           2           3
           2           4           5
           3           5           6
           2           3           5
           ];
ok = any(tri(:)==tri0(:));
