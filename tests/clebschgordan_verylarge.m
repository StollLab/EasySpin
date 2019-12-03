function [err,data] = test(opt,olddata)

% very large J values
%======================================================
a(1) = clebschgordan([5000 25],[4000 -15],[1000 10]);
b(1) = -0.0621923408158706562064638785093;
err = ~areequal(a,b,1e-12,'abs');

data = [];
