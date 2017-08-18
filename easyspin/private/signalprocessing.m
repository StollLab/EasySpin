function [ ProcessedSignal ] = signalprocessing(TimeAxis,RawSignal,Opt)
nDetectionOperators = length(Opt.DetectionOperators);

DownConversionFrequencies = zeros(1,nDetectionOperators);

try
if isfield(Opt,'DownConversionFrequency') &&  ~isempty(Opt.DownConversionFrequency)
  nDownConversionFrequencies = length(Opt.DownConversionFrequency);
  if nDownConversionFrequencies == 1
    DownConversionFrequencies(1:nDetectionOperators) = Opt.DownConversionFrequency;
  else
    DownConversionFrequencies(1:nDownConversionFrequencies) = Opt.DownConversionFrequency;
  end
end

if iscell(RawSignal)
  ProcessedSignal = cell(size(RawSignal));
  nPoints = length(RawSignal);
else
  ProcessedSignal = zeros(size(RawSignal));
  nPoints = size(RawSignal,1);
  DCTimeAxis = TimeAxis;
end

for iPoint = 1 : nPoints
  if iscell(RawSignal)
    Traces = zeros(size(RawSignal{iPoint}));
    DCTimeAxis = TimeAxis{iPoint};
  end
  for iTrace = 1 : nDetectionOperators
    if iscell(RawSignal)
      RfSignal = RawSignal{iPoint}(iTrace,:);
    else
      RfSignal = squeeze(RawSignal(iPoint,iTrace,:));
    end
    
    if max(imag(RfSignal)) < 1e-10
      RfSignal = real(RfSignal);
    elseif max(real(RfSignal)) < 1e-10
      RfSignal = imag(RfSignal);
    else
      RfSignal = 2*RfSignal;
    end
    
    if DownConversionFrequencies(iTrace) ~= 0
      
      if isreal(RfSignal)
        [~, DCSignal] = rfmixer(DCTimeAxis,RfSignal,DownConversionFrequencies(iTrace),'IQdemod');
        if iscell(RawSignal)
          Traces(iTrace,:) = real(DCSignal);
        else
          ProcessedSignal(iPoint,iTrace,:) = real(DCSignal);
        end
      else
        [~, DCSignal] = rfmixer(DCTimeAxis,RfSignal,DownConversionFrequencies(iTrace),'IQshift');
        if iscell(RawSignal)
          Traces(iTrace,:) = DCSignal;
        else
          ProcessedSignal(iPoint,iTrace,:) = DCSignal;
        end
      end
      
    else
      if iscell(RawSignal)
        Traces(iTrace,:) = RfSignal;
      else
        ProcessedSignal(iPoint,iTrace,:) = RfSignal;
      end
      
    end
  end
  if iscell(RawSignal)
    ProcessedSignal{iPoint} = Traces;
  end
end

catch
  ProcessedSignal = RawSignal;
  
  warning('The down conversion of the signal was not successful and created an error. The raw signal in the simulation frame was returned. Please try down conversion manually.')
end

end

