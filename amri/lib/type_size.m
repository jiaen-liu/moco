function bytesperword=type_size(type)

%+
% NAME:
%     
%#
% FUNCTION:
%     Return the size of a word in bytes. Tyep can be an IDL type
%     code or a matlab class.
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
%     #1    Type 'logical' has a size of 1 byte per word.
%           JAdZ, January 31, 2012
%-

% return bytes per word for this type of data
switch type
	case 'uint8',
		bytesperword=1;
	case 'int16',
		bytesperword=2;
	case 'int32',
		bytesperword=4;
	case 'single',
		bytesperword=4;
	case 'double',
		bytesperword=8;
	case 'char',
		bytesperword=1;
	case 'struct',
		bytesperword=0;
	case 'uint16',
		bytesperword=2;
	case 'uint32',
		bytesperword=4;
	case 'int64',
		bytesperword=8;
	case 'uint64',
		bytesperword=8;
	case 'int8',
		bytesperword=1;
	case 'logical',
		bytesperword=1;
	case 'cell',
		bytesperword=0;
	case 'function_handle',
		bytesperword=0;
	case 'numeric',
		bytesperword=0;
	case 1,
		bytesperword=1;
	case 2,
		bytesperword=2;
	case 3,
		bytesperword=4;
	case 4,
		bytesperword=4;
	case 5,
		bytesperword=8;
	case 6,
		bytesperword=8;
	case 7,
		bytesperword=1;
	case 8,
		bytesperword=0;
	case 9,
		bytesperword=16;
	case 10,
		bytesperword=0;
	case 11,
		bytesperword=0;
	case 12,
		bytesperword=2;
	case 13,
		bytesperword=4;
	case 14,
		bytesperword=8;
	case 15,
		bytesperword=8;
	otherwise,
		bytesperword=0;
end
