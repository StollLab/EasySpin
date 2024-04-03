% Plots a CW EPR spectrum (or spectra) simulated using pepper/chili/garlic
function cwepr_plot(xAxis,spec,Exp)

if Exp.FrequencySweep
  if xAxis(end)<1
    xAxis = xAxis*1e3;
    xl = 'frequency (MHz)';
  else
    xl = 'frequency (GHz)';
  end
  titlestr = sprintf('%0.8g mT',Exp.Field);
else
  if xAxis(end)<10000
    xl = 'magnetic field (mT)';
  else
    xAxis = xAxis*1e3;
    xl = 'magnetic field (T)';
  end
  titlestr = sprintf('%0.8g GHz',Exp.mwFreq);
end

crystalRotation = isfield(Exp,'SampleRotation') && ...
  ~isempty(Exp.SampleRotation) && numel(Exp.SampleRotation{2})>1;

if crystalRotation && (size(spec,1)==numel(Exp.SampleRotation{2}))
  h = stackplot(xAxis,spec,'maxabs',2,'');
  set(h,'Color',h(1).Color);
else
  plot(xAxis,spec);
end
axis tight
xlabel(xl);
ylabel('intensity (arb.u.)');
title(titlestr);

end
