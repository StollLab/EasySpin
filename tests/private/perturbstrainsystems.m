function [S,c] = perturbstrainsystems()

% Spin systems with g, A and D strains for the perturbation-theory strain
% tests resfields_perturb_strain_fd and resfreqs_perturb_strain_fd

ang = [0.3 0.7 1.1];
c = 0;
% g principal values and frame, no nuclei
c=c+1; S{c} = struct('g',[2 2.1 2.2],'gFrame',ang,'StrainPars',{{'g(1)','g(3)','gFrame(2)'}},'StrainFWHM',[0.01 0.02 0.05]);
% correlated g and A strains, A frame, HStrain
c=c+1; S{c} = struct('Nucs','63Cu','g',[2.05 2.1 2.25],'gFrame',ang,'A',[60 80 500],'AFrame',[0.5 0.2 0.1], ...
  'StrainPars',{{'g(3)','A(3)','AFrame(3)'}},'StrainFWHM',[0.02 30 0.1],'StrainCorr',[-0.8 0 0],'HStrain',[10 10 20]);
% axial, isotropic and symmetric A
c=c+1; S{c} = struct('Nucs','1H,1H','A',[5 8; 3 9],'StrainPars',{{'A(3)','A(2)'}},'StrainFWHM',[2 1]);
c=c+1; S{c} = struct('Nucs','1H,14N','A',[10 30],'StrainPars',{{'A(2)'}},'StrainFWHM',5);
c=c+1; S{c} = struct('Nucs','1H','g',[2 2.1 2.2],'A',[5 8 12 1 2 3],'StrainPars',{{'g(2)','A(5)'}},'StrainFWHM',[0.02 2]);
% g strain with non-symmetric full A, A strain with non-symmetric full g
c=c+1; S{c} = struct('Nucs','1H','g',[2 2.1 2.2],'gFrame',ang,'A',[10 2 1; -1 15 3; 0.5 -2 20], ...
  'StrainPars',{{'g(1)','gFrame(1)'}},'StrainFWHM',[0.02 0.1]);
c=c+1; S{c} = struct('Nucs','1H','g',[2 0.02 0.01; -0.01 2.1 0.03; 0.02 0 2.2],'A',[10 15 20],'AFrame',ang, ...
  'StrainPars',{{'A(2)','AFrame(2)'}},'StrainFWHM',[3 0.1]);
% S = 1: [D E] and frame, principal values, symmetric D, g strain with full D
c=c+1; S{c} = struct('S',1,'D',[1000 200],'DFrame',ang,'StrainPars',{{'D(1)','D(2)','DFrame(1)'}},'StrainFWHM',[50 20 0.05]);
c=c+1; S{c} = struct('S',1,'D',[-300 -100 500],'StrainPars',{{'D(1)','D(3)'}},'StrainFWHM',[30 40]);
c=c+1; S{c} = struct('S',1,'D',[-300 -100 400 50 30 20],'StrainPars',{{'D(4)','D(1)'}},'StrainFWHM',[20 30]);
c=c+1; S{c} = struct('S',1,'g',[2 2.1 2.2],'D',[100 20 10; 20 -50 5; 10 5 300],'StrainPars',{{'g(2)','g(3)'}},'StrainFWHM',[0.01 0.02]);
% S = 3/2 (central transition), D strain with nucleus
c=c+1; S{c} = struct('S',3/2,'D',[800 100],'StrainPars',{{'D(1)'}},'StrainFWHM',100);
c=c+1; S{c} = struct('S',3/2,'D',[500 50],'DFrame',ang,'Nucs','1H','A',[20 30 40],'StrainPars',{{'D(1)','A(2)'}},'StrainFWHM',[60 5]);
% strain modes
c=c+1; S{c} = struct('Nucs','1H','g',[2 2.1 2.2],'A',[10 20 30],'StrainPars',{{'g(1)','A(3)'}},'StrainModes',[0.01 5; 0 3]);
% all strain FWHM zero (including Q), with HStrain
c=c+1; S{c} = struct('Nucs','14N','g',2,'A',[10 20 30],'Q',1,'StrainPars',{{'g','Q'}},'StrainFWHM',[0 0],'HStrain',[5 5 8]);
