function ok = test()

% Without frequency translation, signals are returned unchanged, also if
% they are purely imaginary or purely real

t = 0:0.001:0.2;  % µs
Freq = 150;  % MHz

ImagSignal = 1i*exp(-t/0.05);
RealSignal = exp(-t/0.05);
ComplexSignal = exp(2i*pi*Freq*t - t/0.05);

% one acquisition point, three detection operators
RawSignal = reshape([ImagSignal; RealSignal; ComplexSignal],[1 3 numel(t)]);
Signal = runprivate('signalprocessing',t,RawSignal,[0 0 0]);
ok(1) = areequal(Signal(:,1),ImagSignal(:),1e-12,'abs');
ok(2) = areequal(Signal(:,2),RealSignal(:),1e-12,'abs');
ok(3) = areequal(Signal(:,3),ComplexSignal(:),1e-12,'abs');

% With frequency translation, a purely imaginary signal stays imaginary
RawSignal = 1i*cos(2*pi*Freq*t).';
Signal = runprivate('signalprocessing',t,RawSignal,-Freq/1e3);
ok(4) = max(abs(real(Signal))) < 1e-12 && max(abs(imag(Signal))) > 0.1;
