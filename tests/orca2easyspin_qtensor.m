function [err,data] = test(opt,olddata)

% Test whether Q tensor is read in correctly from ORCA property file
%--------------------------------------------------------------------
Sys = orca2easyspin('orca/hydroxyl_098.oof');
Q0 = [-0.0644363 -0.0896791 0.1541154]; % MHz
err = any(abs(Sys.Q-Q0)>1e-5);
data = [];
