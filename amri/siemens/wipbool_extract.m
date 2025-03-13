function wipb=wipbool_extract(tfree)

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
%     July 12, 2022
%     Jacco de Zwart
%     LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%     E-mail: Jacco.deZwart@nih.gov
%
%; INFORMATION:
%     Character map used to encode 64 bits, aka 6 bools:
%         045-045: -        21
%         048-057: 0..9     22..31
%         065-090: A..Z     32..57
%         095-095: _        58
%         097-122: a..z     59..63,0..20 --> f=0 ---> bits=000000, for backward compatibility
%
% MODIFICATION HISTORY:
%    2022-07-14, JAdZ
%        Removed residual debugging print statements.
%-

% definitions
dense_tfree_offset=97;

% create a temporary of the correct length - 6 bytes needed for each char > dense_tfree_offset
tsiz=size_of(tfree);
if tsiz > dense_tfree_offset
	tsiz=tsiz+(size_of(tfree)-dense_tfree_offset)*5;
end
wipb=zeros(1,tsiz,'uint8');
wipb(1:size_of(tfree))=uint8(tfree == 't');
if tsiz > dense_tfree_offset
	for i=(dense_tfree_offset+1):size_of(tfree)

		% covert to 0..63 range
		cval=tfree(i);
		done=0;
		if (cval == 45)
			cval=0;
			done=1;
		end
		if (cval >= 48) && (cval <= 57)
			cval=cval-47;
			done=1;
		end
		if (cval >= 65) && (cval <= 90)
			cval=cval-54;
			done=1;
		end
		if (cval == 95)
			cval=cval-58;
			done=1;
		end
		if (cval >= 97) && (cval <= 122)
			cval=cval-59;
			done=1;
		end
		if (done == 0)
			error("wipbool_extract()> ERROR: Unexpected characters in tFree string!");
		end

		% shift by 21 and mod so that 'f' is 0
		cval=cval+21;
		cval=mod(cval,64);

		% convert to 6 bits
		bits=dec2bin(cval);

		% store in the result array
		while (size_of(bits) < 6)
			bits=strcat('0',bits);
		end
		offs=6*(i-dense_tfree_offset-1)+dense_tfree_offset+1;
		for j=0:5
			wipb(offs+j)=int8(str2num(bits(6-j)));
		end
	end
end
