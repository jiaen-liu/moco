function type=type_of(data)

%+
% NAME:
%     Return the IDL data type.
%#
% FUNCTION:
%     
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
%     #1    Can now properly handle complex data.
%           JAdZ, April 2, 2008
%     #2    Defined both 'logical' and 'int8' as byte datatype (type=1).
%           Cell-type data is tested for being a 'string array'.
%           JAdZ, January 31, 2012
%     #3    Added type 'string'.
%           JAdZ, September 24, 2018
%-

mtltype=class(data);
switch mtltype
	case 'uint8',
		type=1;
	case 'int8',
		type=1;
	case 'logical',
		type=1;
	case 'int16',
		type=2;
	case 'int32',
		type=3;
	case 'single',
		if isreal(data),
			type=4;
		else
			type=6;
		end
	case 'double',
		if isreal(data),
			type=5;
		else
			type=9;
		end
	case 'char',
		type=7;
	case 'string',
		type=7;
	case 'struct',
		type=8;
	case 'uint16',
		type=12;
	case 'uint32',
		type=13;
	case 'int64',
		type=14;
	case 'uint64',
		type=15;
	case 'cell',
		if iscellstr(data)
			type=7;
		else
			type=0;
		end
	case 'function_handle',
		type=0;
	case 'numeric',
		type=0;
	otherwise,
		type=0;
end
