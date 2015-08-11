clear;

%profile on;

%--------------------------------------------------------------------------
% call my SLE program given the system and experimental structures:
%    Sys:
%        S -- total spin
%        g -- g-tensor
%        A -- hyperfine tensor(s)
%        D -- ZFS tensor
%    Exp:
%        mw   -- microwave frequency
%        B    -- static magnetic field
%        temp -- temperature
%--------------------------------------------------------------------------

%sle(Sys,Exp);

Spins = 3;
gTensor = 3;
hfTensor = 5;
zfsTensor = 3;

switch Spins
  case 1 % free electron
    Sys.S = 1/2;
  case 2 % two uncoupled electrons
    Sys.S = [1/2 1/2];
  case 3 % triplet
    Sys.S = 1;
  case 4 % gadolinium
    Sys.S = 7/2;
  case 5 % manganese
    Sys.S = 5/2;
end

switch gTensor
  case 1 % organic radical (nitronyl nitroxide)
    %Sys.g = diag([2.01216 2.00759 2.00336]);
    Sys.g = [2.01216 2.00759 2.00336];
  case 2 % highly anisotropic
    Sys.g = [2.05 2.025 2.0];
  case 3 % isotropic
    giso = 2;%1.997;
    Sys.g = giso;
    Sys.g = [giso giso giso];
  case 4 % axial
    gpar = 2.08;
    gperp = 2.06;
    Sys.g = diag([gperp gperp gpar]);
end

%Sys.A = cell(numel(find(SpinType)),1)
switch hfTensor
  case 1 % nitroxide
    Sys.Nucs = '14N';
    Sys.A = [20 30 50];
  case 2 % nitronyl nitroxide
    Sys.Nucs = '14N,14N';
    Sys.A = [10 10 40; 10 10 40];
  case 3 % manganese
    Sys.Nucs = '15N';
    %Sys.A = [50 50 90];
    Sys.A = 10;
  case 4 % four nuclei
    Sys.Nucs = '14N,14N,1H,1H';
    Sys.A = [10 10 40; 10 10 40; 5 5 10; 5 5 10];
  case 5 % none
end

switch zfsTensor
  case 1 % large
    Sys.D = 500; 
  case 2 % medium
    Sys.D = 250;
  case 3 % small
    Sys.D = 50;
  case 4 % none
end

Sys.tcorr = 100e-9;
%Sys.lwpp = 2;

%Exp.mwFreq = 9.5;
Exp.Field = 334;
Exp.Range = [9.28 9.42];
%Exp.nPoints = 2048;
%Exp.CenterSweep = [Exp.Field*bmagn/(planck/(2*pi)), 20];

Opt.LiouvMethod = 'general';
Opt.LLKM = [24 7 6 6];
%[Sys,Exp] = js2es(Sys,Exp);

%sle(Sys,Exp);



return