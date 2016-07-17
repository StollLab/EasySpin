% An partially oriented Cu(II) complex
%==========================================================================
% The molecule is assumed to be planar, with the z axis
% of the molecular reference frame perpendicular to the plane, just
% as in, e.g., CuTPP.

clear, clf

% Paramagnetic complex
Sys = struct('g',[2 2 2.2],'lwpp',0.8);
Sys = nucspinadd(Sys,'63Cu',[50 50 400]);

% Experimental parameters
Exp = struct('mwFreq',9.5,'Range',[280 350]);

% (1) Isotropic orientational distribution
[B,spec1] = pepper(Sys,Exp);
% (2) Partially oriented, B0 approx. in molecular xy plane
Exp.Ordering = -2;
[B,spec2] = pepper(Sys,Exp);
% (3) Strongly oriented, B0 along molecular z axis
Exp.Ordering = 4;
[B,spec3] = pepper(Sys,Exp);

plot(B,spec1,'k',B,spec2,'g',B,spec3,'r')
legend('isotropic distribution','preferential in xy plane','preferential along z axis');
