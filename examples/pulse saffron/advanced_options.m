% Advanced options for spin dynamics (saffron)
%==========================================================================
% demonstration of some of the advanced options in saffron for the spin
% dynamics of a single orientation and how to use them

clear

% Default Spin System
Sys.S = 1/2;
Sys.T1 = 1; % µs
Sys.T2 = 0.5; % µs

% Pulse Definitions
Rectangular.Type = 'rectangular';
Rectangular.tp = 0.02; % µs
Rectangular.Flip = pi/2; %rad

Adiabatic.Type = 'quartersin/linear';
Adiabatic.Flip = pi/2; % rad
Adiabatic.tp = 0.2; % µs
Adiabatic.trise = 0.05; % µs
Adiabatic.Qcrit = 5; % critical adiabaticity
Adiabatic.Frequency = [-50 50]; % frequency band,MHz

% A default Experiment/Sequence
Exp.mwFreq = 33.5; % GHz
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree at 33.5 GHz
Exp.MolFrame = [0 0 0]; % single orientation, molecular frame = lab frame
Exp.DetWindow = 'sequence'; % detect during the whole sequence
Exp.DetOperator = {'z1'};

%% 1) Preparing a state and setting the equilibrium state
% An ideal pi pulse (no pulse length) at the start inverts the thermal
% state -Sz to +Sz. During the free evolution, the spin relaxes to the
% equilibrium state Sys.eqState. Relaxation is switched off during the
% pulses (one entry of Opt.Relaxation per element of Exp.Sequence).

Sys_ = Sys;
Sys_.eqState = -sop(Sys_.S,'z');

IdealPi.Flip = pi; % rad, ideal pulse

Exp_ = Exp;
Exp_.Sequence = {IdealPi Rectangular 5}; % ideal pi pulse, pulse and 5 µs free evolution

Opt_.Relaxation = [0 0 1];

[TimeAxis, Signal] = saffron(Sys_,Exp_,Opt_);

% plotting
figure(1)
clf
plot(TimeAxis*1000,real(Signal))
xlabel('t (ns)')
ylabel('<S_i>')
legend(Exp.DetOperator)

%% 2) circular mwPolarization
% If mwPolarization is 'circular', the excitation operator takes the form:
% real(IQ)*Sx + imag(IQ)*Sy
% In a shifted simulation frame (the default), pulses keep only the field
% component that rotates with the spins, so linear and circular
% polarization give the same result. They differ only in the lab frame,
% where the counter-rotating component of a linearly polarized field
% causes the Bloch-Siegert shift and oscillations at twice the microwave
% frequency. These effects scale with the ratio of the pulse amplitude to
% the microwave frequency, so this example uses a strong pulse at a low
% frequency and simulates in the lab frame (Opt.SimFreq = 0).

Sys_ = [];
Sys_.S = 1/2;

Strong.Type = 'rectangular';
Strong.tp = 0.005; % µs
Strong.Flip = pi; % rad

Exp_ = [];
Exp_.mwFreq = 0.3; % GHz
Exp_.Field = mhz2mt(300);  % mT, resonance field for g = gfree at 0.3 GHz
Exp_.MolFrame = [0 0 0];
Exp_.Sequence = {Strong};
Exp_.DetWindow = 'sequence';
Exp_.DetOperator = {'z1'};

Opt_ = [];
Opt_.SimFreq = 0;  % lab frame

[TimeAxis, Regular] = saffron(Sys_,Exp_,Opt_);

Exp_.mwPolarization = 'circular';

[~, ComplexExOp] = saffron(Sys_,Exp_,Opt_);

% plotting
figure(2)
clf
hold on
plot(TimeAxis*1000,Regular)
plot(TimeAxis*1000,ComplexExOp)
xlabel('t (ns)')
ylabel('<S_i>')
legend('linear mwPolarization','circular mwPolarization')

%% 3) Bloch vector during the sequence
% Detecting Sx, Sy and Sz during the whole sequence follows the Bloch
% vector of the spin. Sx and Sy are given in the frame rotating at
% Exp.DetFreq (default Exp.mwFreq). The adiabatic pulse inverts the spin
% (with Qcrit given, Flip is ignored), and relaxation (Sys.T1, Sys.T2)
% shrinks the Bloch vector.

Exp_ = Exp;
Exp_.Sequence = {Adiabatic 0.5};
Exp_.DetOperator = {'x1','y1','z1'};

[TimeAxis, Bloch] = saffron(Sys,Exp_);

% plotting
figure(3)
clf
plot3(real(Bloch(:,1)),real(Bloch(:,2)),real(Bloch(:,3)));
xlabel('<S_x>')
ylabel('<S_y>')
zlabel('<S_z>')
axis equal
grid on
