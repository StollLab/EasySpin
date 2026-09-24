% xml2struct       Convert XML file into a MATLAB structure
%
%     s = xml2struct(fileName)
%
%  Input:
%     fileName    file name of XML file (extension .xml can be omitted)
%
%  Output:
%     s           A hierarchical MATLAB structure representing the XML tree
%
%  A file containing:
%    <XMLname attrib1="Some value">
%      <Element>Some text</Element>
%      <Data attrib2="2">Some more text</Data>
%      <Data attrib3="8.3" attrib4="1">Even more text</Data>
%    </XMLname>
%
%  will produce:
%     s.XMLname.Attributes.attrib1 = 'Some value';
%     s.XMLname.Element.Text = 'Some text';
%     s.XMLname.Data{1}.Attributes.attrib2 = '2';
%     s.XMLname.Data{1}.Text = 'Some more text';
%     s.XMLname.Data{2}.Attributes.attrib3 = '8.3';
%     s.XMLname.Data{2}.Attributes.attrib4 = '1';
%     s.XMLname.Data{2}.Text = 'Even more text';
%
%  All values are char row vectors. An element that occurs once is a struct,
%  an element that occurs several times is a cell array of structs.
%
%  Text, CDATA, and comment content of an element is stored in the fields
%  Text, CDATA, and Comment, respectively, alongside any child elements.
%  Whitespace-only content is omitted. Processing instructions are skipped.
%
%  In element and attribute names, the following characters are substituted:
%   '-'     '_dash_'
%   ':'     '_colon_'
%   '.'     '_dot_'
%  Any remaining characters not allowed in MATLAB field names are replaced
%  by '_', and names are truncated to namelengthmax.
%
%  The file is parsed with the MATLAB DOM API (matlab.io.xml.dom, R2021a and
%  later), which does not require Java. Files with a DOCTYPE declaration are
%  rejected, as a protection against XXE attacks.

% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012
% Modified by Stefan Stoll, University of Washington, Feb 2018, Sep 2026

function s = xml2struct(fileName)

if nargin < 1
  help(mfilename);
  return
end

fileName = char(fileName);

% Check for existence, add .xml extension if omitted
if ~isfile(fileName) && ~endsWith(fileName,'.xml','IgnoreCase',true)
  fileName = [fileName '.xml'];
end
if ~isfile(fileName)
  error('The file %s could not be found.',fileName);
end

% Read the xml file
xDoc = parseFile(matlab.io.xml.dom.Parser,fileName);

% Parse xDoc into a MATLAB structure
s = parseChildNodes(xDoc);

end

% ----- Subfunction parseChildNodes -----
function [children,text] = parseChildNodes(theNode)
% Recurse over node children. Returns child elements in children, and
% text, CDATA, and comment content in text.

children = struct;
text = struct;
if ~hasChildNodes(theNode)
  return
end

childNodes = getChildNodes(theNode);
for count = 1:getLength(childNodes)
  theChild = item(childNodes,count-1);

  if isa(theChild,'matlab.io.xml.dom.Element')
    name = validName(getNodeName(theChild));
    element = parseElement(theChild);
    % XML allows the same element to occur multiple times,
    % put each occurrence in a different cell
    if isfield(children,name)
      if ~iscell(children.(name))
        children.(name) = {children.(name)};
      end
      children.(name){end+1} = element;
    else
      children.(name) = element;
    end
    continue
  end

  % CDATASection is a subclass of Text, so it has to be checked first
  if isa(theChild,'matlab.io.xml.dom.CDATASection')
    textfield = 'CDATA';
  elseif isa(theChild,'matlab.io.xml.dom.Text')
    textfield = 'Text';
  elseif isa(theChild,'matlab.io.xml.dom.Comment')
    textfield = 'Comment';
  else
    continue % skip processing instructions etc.
  end

  str = rowchar(getTextContent(theChild));
  if all(isspace(str)), continue; end
  if isfield(text,textfield)
    % <element>Text <!--Comment--> More text</element>: append the text
    text.(textfield) = [text.(textfield) str];
  else
    text.(textfield) = str;
  end

end

end

% ----- Subfunction parseElement -----
function element = parseElement(theNode)
% Create structure of element: child elements, text, and attributes.

[element,text] = parseChildNodes(theNode);

% Add text content alongside child elements
textfields = fieldnames(text);
for k = 1:numel(textfields)
  element.(textfields{k}) = text.(textfields{k});
end

% Store the (possibly whitespace-only) text of elements without other content
if isempty(fieldnames(element))
  element.Text = rowchar(getTextContent(theNode));
end

attributes = parseAttributes(theNode);
if ~isempty(fieldnames(attributes))
  element.Attributes = attributes;
end

end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = struct;
if hasAttributes(theNode)
  theAttributes = getAttributes(theNode);
  for count = 1:getLength(theAttributes)
    attrib = item(theAttributes,count-1);
    attributes.(validName(getName(attrib))) = rowchar(getValue(attrib));
  end
end

end

% ----- Subfunction validName -----
function name = validName(name)
% Convert XML name to valid MATLAB field name.

name = char(name);
name = strrep(name, '-', '_dash_');
name = strrep(name, ':', '_colon_');
name = strrep(name, '.', '_dot_');
name = matlab.lang.makeValidName(name);
name = name(1:min(end,namelengthmax));

end

% ----- Subfunction rowchar -----
function str = rowchar(str)
% Convert to char row vector (also for empty strings).
str = reshape(char(str),1,[]);
end
