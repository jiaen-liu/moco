function res=swap_endian(data)

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
%     July 21, 2004
%     Jacco de Zwart
%     LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%     E-mail: Jacco.deZwart@nih.gov
%
% MODIFICATION HISTORY:
%     #1    Now uses typecast in stead of cast.
%           JAdZ, December 16, 2005
%-

% reform data to single-dimensional array
elms=n_elements(data);
res=data;
class(data);
wordsi=type_size(class(data));
if [ wordsi == 0 ]
	disp('Unable to byte-swap this data type');
	return
end

% loop over the number of elements and bytes
newword=uint8(zeros(1,wordsi));
for i=1:elms,
	oldword=typecast(data(i),'uint8');
	for j=1:wordsi,
		newword(wordsi-j+1)=oldword(j);
	end
	res(i)=typecast(newword,class(data));
end

% convert res and data back to original dimensions
