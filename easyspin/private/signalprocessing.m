function [ ProcessedSignal ] = signalprocessing(TimeAxis,RawSignal,Opt)
% Number of detection operators
nDetectionOperators = length(Opt.DetectionOperators);

% Setup vector with down conversion frequencies
TranslationFrequencies = zeros(1,nDetectionOperators);

% Trys down conversion, if it fails for any reason, the raw signal is
% returned. Errors are usually when the down conversion frequency is very
% wrong, or a non oscillating signal is to be down converted (Sz)
try
  % This transcribes the down conversion frequencies into the corresponding
  % vector. If only one down conversion frequency is provided, this is
  % being used for all detection operators. 
  %%%%%%%%%%%%%%%%%%
  % We might want to remove this feature though, too much of an assumption
  if isfield(Opt,'FreqTranslation') &&  ~isempty(Opt.FreqTranslation)
    nDownConversionFrequencies = length(Opt.FreqTranslation);
    if nDownConversionFrequencies == 1
      TranslationFrequencies(1:nDetectionOperators) = Opt.FreqTranslation;
    else
      TranslationFrequencies(1:nDownConversionFrequencies) = Opt.FreqTranslation;
    end
  end
  
  % Recognizing the type of the input, wheter it is a cell or numeri array
  % and gets the number of acquisition points
  % For the cell array each element has to processed individually
  if iscell(RawSignal)
    ProcessedSignal = cell(size(RawSignal));
    nPoints = length(RawSignal);
  else
    ProcessedSignal = zeros(size(RawSignal));
    nPoints = size(RawSignal,1);
    DCTimeAxis = TimeAxis;
  end
  
  % loop over all acquired data points
  for iPoint = 1 : nPoints
    if iscell(RawSignal)
      % Gets the size of the traces for the current data point and the time
      % axis if the RawSignal is a CellArray
      Traces = zeros(size(RawSignal{iPoint}));
      DCTimeAxis = TimeAxis{iPoint};
    end
    % Loop over the Number of Detection Operators
    for iTrace = 1 : nDetectionOperators
      % Load RawSignal of the current Detection Operator
      if iscell(RawSignal)
        RfSignal = RawSignal{iPoint}(iTrace,:);
      else
        RfSignal = squeeze(RawSignal(iPoint,iTrace,:));
      end
      
      % Ensures that signal is purely real or imaginary, if it is supposed 
      % to be purely real/imaginary
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Compare to the magnitude
      
      if max(imag(RfSignal)) < 1e-10
        RfSignal = real(RfSignal);
      elseif max(real(RfSignal)) < 1e-10
        RfSignal = imag(RfSignal);
      else
        % If signal is complex (S+ or S-), the normalization is wrong by a
        % factor of two and hence needs normalization
        RfSignal = 2*RfSignal;
      end
      
      % Does the downconversion, if no down conversion is requested for the
      % current detection operator, the cleaned up signal is written to
      % ProcessedSignal
      if TranslationFrequencies(iTrace) ~= 0
        % Depending on the type of trace, different down conversion types
        % are required (see rfmixer for details)
        if isreal(RfSignal)
          % Mixing if signal is real
          [~, DCSignal] = rfmixer(DCTimeAxis,RfSignal,TranslationFrequencies(iTrace),'IQdemod');
          if iscell(RawSignal)
            Traces(iTrace,:) = real(DCSignal);
          else
            ProcessedSignal(iPoint,iTrace,:) = real(DCSignal);
          end
        else
          % Mixing if signal is complex
          [~, DCSignal] = rfmixer(DCTimeAxis,RfSignal,TranslationFrequencies(iTrace),'IQshift');
          if iscell(RawSignal)
            Traces(iTrace,:) = DCSignal;
          else
            ProcessedSignal(iPoint,iTrace,:) = DCSignal;
          end
        end
        
      else
        % stores original cleaned up signal
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
  % If something goes wrong during down conversion, the RawSignal is
  % returned
  ProcessedSignal = RawSignal;
  
  warning('The down conversion of the signal was not successful and created an error. The raw signal in the simulation frame was returned. Please try down conversion manually.')
end

end

