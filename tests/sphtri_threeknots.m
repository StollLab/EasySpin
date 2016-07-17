function [err,data] = test(opt,olddata)

tri = sphtri('D2h',3);
tri0 = [           
           1           2           3
           2           4           5
           3           5           6
           2           3           5
           ];
err = any(tri(:)~=tri0(:));
data = [];
