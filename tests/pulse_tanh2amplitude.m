function ok = test(opt)

% Tanh pulse amplitude shape
tp = 0.200;  % pulse length, µs
tr = 0.020;  % rise time, µs

Params.tp = tp;
Params.Type = 'tanh2';
Params.trise = tr;
Params.Amplitude = 1;

[t,IQ] = pulse(Params);

Amplitude = coth(tp/2/tr)^4 * tanh(t/tr).^2 .* tanh((tp-t)/tr).^2;
IQref = Amplitude;

if opt.Display
  plot(t,IQref,t,IQ);
  legend('reference','pulse()');
  xlabel('t (µs)');
  ylabel('amplitude');
end

ok = areequal(IQ,IQref,1e-12,'abs');
