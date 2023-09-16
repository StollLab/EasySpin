% Plots a CW EPR spectrum (or spectra) simulated using pepper/chili/garlic
function cwepr_plot(xAxis,spec,Exp)

cla
if Exp.FrequencySweep
  if xAxis(end)<1
    plot(xAxis*1e3,spec);
    xlabel('frequency (MHz)');
  else
    plot(xAxis,spec);
    xlabel('frequency (GHz)');
  end
  title(sprintf('%0.8g mT',Exp.Field));
else
  if xAxis(end)<10000
    plot(xAxis,spec);
    xlabel('magnetic field (mT)');
  else
    plot(xAxis/1e3,spec);
    xlabel('magnetic field (T)');
  end
  title(sprintf('%0.8g GHz',Exp.mwFreq));
end
axis tight
ylabel('intensity (arb.u.)');


end

