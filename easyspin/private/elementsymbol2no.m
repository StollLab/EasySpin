% Convert element symbol to element number.
% Example: 'Cu' or "Cu" -> 29

function no = elementsymbol2no(symbol)

% Convert strings to char arrays
if isstring(symbol)
  symbol = char(symbol);
end

% Crop and pad symbol
if numel(symbol)==1
  symbol(2) = ' ';
elseif numel(symbol)>2
  symbol  = symbol(1:2);
end

ElementSymbols = [...
  'H He'...
  'LiBeB C N O F Ne'...
  'NaMgAlSiP S ClAr'...
  'K CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKr'...
  'RbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI Xe'...
  'CsBaLaCePrNdPmSmEuGdTbDyHoerTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn'...
  'FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDsRgCn'];

idx = strfind(ElementSymbols,symbol);
no = (idx+1)/2;
