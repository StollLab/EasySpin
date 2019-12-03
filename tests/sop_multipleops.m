function [err,data] = test(opt,olddata)

% Multiple spin operators, syntax check

S = 1;
[Sx,Sy,Sz,Sp] = sop(S,'x','y','z','+');

Sx_ = sop(S,'x');
Sy_ = sop(S,'y');
Sz_ = sop(S,'z');
Sp_ = sop(S,'+');

thr = 1e-12;
err = ~areequal(Sx,Sx_,thr,'abs') || ~areequal(Sy,Sy_,thr,'abs') || ...
  ~areequal(Sz,Sz_,thr,'abs') || ~areequal(Sp,Sp_,thr,'abs');

data = [];
