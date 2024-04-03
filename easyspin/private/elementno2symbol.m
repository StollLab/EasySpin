% Convert element number to element symbol.
% Example: 29 -> 'Cu'

function symbol = elementno2symbol(no)

if no<1 || no>112
  symbol = '';
  return
end

ElementSymbols = [...
  'H He'...
  'LiBeB C N O F Ne'...
  'NaMgAlSiP S ClAr'...
  'K CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKr'...
  'RbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI Xe'...
  'CsBaLaCePrNdPmSmEuGdTbDyHoerTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn'...
  'FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDsRgCn'];
idx = no*2 - 1;

symbol = ElementSymbols(idx:idx+1);

% Remove trailing space for H, C, N, etc.
if symbol(2)==' '
  symbol = symbol(1);
end
