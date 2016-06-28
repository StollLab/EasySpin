% First order sech/tanh and higher order sech inversion pulses
%==========================================================================
% In this example, offset-independent adiabaticity inversion pulses with
% amplitude modulation functions determined by sech of different order are
% computed and their amplitude and frequency modulation functions as well
% as excitation profiles are compared.

clear, clf

nmax = 8;

% Loop over sech pulses with different order n
for n = 1:nmax
  
  % Pulse parameter definitions
  %-------------------------------------------------------------
  Par.Type = 'sech/uniformQ'; % pulse shape, sech amplitude modulated pulse
  % with uniform adiabaticity (n = 1 corresponds
  % to the sech/tanh pulse)
  Par.tp = 0.200; % pulse length, µs
  Par.BW = 120; % pulse bandwidth, MHz
  Par.beta = 6; % truncation parameter, used as (beta/tp)
  Par.n = n;
  Par.Flip = pi; % pulse flip angle
  
  [t{n},y{n},p{n},m{n}] = pulse(Par);
  
  % Plot
  %-------------------------------------------------------------  
  cc = jet(nmax);
  subplot(3,1,1)
  hold on; box on;
  title('Amplitude modulation')
  plot(t{n},m{n}.A,'Color',cc(n,:));
  xlabel('t (\mus)')
  ylabel('\nu_1 (MHz)')
  axis tight

  subplot(3,1,2)
  hold on; box on;
  title('Frequency modulation')
  plot(t{n},m{n}.nu,'Color',cc(n,:));
  xlabel('t (\mus)')
  ylabel('\nu (MHz)')
  axis tight
  
  subplot(3,1,3)
  hold on; box on;
  title('Excitation profile')
  if n==1
    line([0 0],[-1 1],'Color',[1 1 1]*0.8);
    line([min(p{n}.offsets) max(p{n}.offsets)],[0 0],'Color',[1 1 1]*0.8);
    line([1 1]*(-Par.BW/2),[-1 1],'Color',[1 1 1]*0.8);
    line([1 1]*(+Par.BW/2),[-1 1],'Color',[1 1 1]*0.8);
  end
  plot(p{n}.offsets,p{n}.Mz,'Color',cc(n,:));
  xlabel('t (\mus)')
  ylabel('M_z/M_0')
  axis tight
  
  lgd{n} = ['n = ',num2str(n)];
  
end
subplot(3,1,1)
legend(lgd)