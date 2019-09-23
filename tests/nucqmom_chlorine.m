function [err,data] = test(opt,olddata)

% Cl values
%======================================================
q = nucqmom('35Cl,37Cl');
%q1 = [-0.08165 -0.06435];
q1 = [-0.0817  -0.0644 ];
err = ~areequal(q,q1,1e-10,'rel');
data = [];
