function ok = test()

% Processing instructions are skipped, not turned into elements

fileName = [tempname '.xml'];
writexml('<?xml version="1.0"?><?xml-stylesheet href="a.xsl"?><a><?pi x?>1</a>',fileName);
cleanup = onCleanup(@()delete(fileName));

s = runprivate('xml2struct',fileName);

ok = isequal(fieldnames(s),{'a'}) && isequal(s.a,struct('Text','1'));

end

%-------------------------------------------------------------------------------
function writexml(xml,fileName)
fid = fopen(fileName,'w');
fwrite(fid,unicode2native(xml,'UTF-8'));
fclose(fid);
end
