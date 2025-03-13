function vers=hr_ideaversion(ulvers)

%+
% NAME:
%     
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
%     September 20, 2018
%     Jacco de Zwart
%     LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%     E-mail: Jacco.deZwart@nih.gov
%
% MODIFICATION HISTORY:
%     2022-07-11, JAdZ
%         Added VE12U.
%     2022-07-15, JAdZ
%         Added a 2nd entry for VE11C (51130001).
%     2022-08-16, JAdZ
%         Merged in variants from Jiaen's version. Now calls error() when the
%         version is unknown.
%-

% Based on the IDL version, the input can be text, hexadecimal or decimal?
%   Check for the hex first, ignore after :, which is added by Jon Polimeni's
%   read_meas_prot.m.
pos=strfind(ulvers,':');
if [numel(pos) > 0]
	vers=extractBetween(ulvers,1,pos(1)-1);
else
	vers=ulvers;
end
vers=string(vers(1));

% analyze hex
switch vers
	case '0xbee332'
		vers='VA25A';
	case '0x1421cf5'
		vers='VB11D';
	case '0x1452a3b'
		vers='VB13A';
	case '0x1483779'
		vers='VB15A';
	case '21710006';
		vers='VB17A';
	case '0x14b44b6'
		vers='VB17A';
	case '0x273bf24'
		vers='VD11D';
	case '0x2765738'
		vers='VD13A';
	case '0x276a554'
		vers='VD13C';
	case '0x276cc66'
		vers='VD13D';
	case '51110009'
		vers='VE11A';
	case '0x30c0783'
		vers='VE11B';
	case '0x30c2e91'
		vers='VE11C';
	case '51150000'
		vers='VE11E';
	case '51180001'
		vers='VE11K';
	case '51130001'
		vers='VE12U';
	case '51280000'
		vers='VE12U';
	otherwise
		error('Siemens IDEA N4 version ''%s'' not yet defined!\n',vers);
end
