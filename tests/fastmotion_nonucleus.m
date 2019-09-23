function [err,data] = test(opt,olddata)

% Linewdiths of hypothetical system without nuclear spins
%======================================================

Sys.g = [2.0 2.1 2.2];
Field = 350;
tcorr = 1e-11;

lw = fastmotion(Sys,Field,tcorr);

lw_correct = 0.41832353;

err = ~areequal(lw,lw_correct,1e-5,'rel');
data = [];
