function class=type2class(idltype,varargin)

%+
% NAME:
%     
%#
% FUNCTION:
%     Return the matlab class to which this idl data type corresponds.
%
% USAGE:
%     
%
% INPUT:
%     
% OUTPUT:
%     
% RETURNS:
%     
%@
% CALLS:
%     
%
% DISCLAIMER AND CONDITIONS FOR USE:
%     Use of this software is at the user's OWN RISK. Functionality
%     is not guaranteed by creator nor modifier(s), if any.
%     This software may be freely copied and distributed. The original 
%     header MUST stay part of the file and modifications MUST be
%     reported in the 'MODIFICATION HISTORY'-field, including the
%     modification date and the name of the modifier.
%
% CREATED:
%     July 20, 2004
%     Jacco de Zwart
%     LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%     E-mail: Jacco.deZwart@nih.gov
%
% MODIFICATION HISTORY:
%     #1    Complex and dcomplex data are now defined as single and double,
%           respectively, and then reformed later.
%           JAdZ, September 27, 2006
%     #2    Keyword fread_raw added, which will return 'type=>type' instead of
%           'type', for use in fread.
%           JAdZ, October 8, 2020
%-

% keywords
fread_raw=0;
p=inputParser;
var=[];
addParameter(p,'fread_raw',fread_raw,@isnumeric);
p.parse(varargin{:});
fread_raw=p.Results.fread_raw;

% set the class type
switch idltype(1)
	case 1;
		class='uint8';
	case 2;
		class='int16';
	case 3;
		class='int32';
	case 4;
		class='single';
	case 5;
		class='double';
	case 6;
		class='single';
	case 7;
		class='char';
	case 8;
		class='struct';
	case 9;
		class='double';
	case 10;
		class='';
	case 11;
		class='';
	case 12;
		class='uint16';
	case 13;
		class='uint32';
	case 14;
		class='int64';
	case 15;
		class='uint64';
	otherwise,
		disp('Unknown data type or class');
		class='';
end

% fread raw mode
if [ fread_raw ~= 0 ]
	class=strcat(class,'=>',class);
end
