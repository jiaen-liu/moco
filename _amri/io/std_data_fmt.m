function dstd=std_data_fmt(data)

%+
% NAME:
%     
%#
% FUNCTION:
%     Converts data to a standard type of byte-type data. See the IDL version
%     of this function for more information.
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
%     January 31 - February 1, 2012
%     Jacco de Zwart
%     LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%     E-mail: Jacco.deZwart@nih.gov
%
% MODIFICATION HISTORY:
%-

% get the IDL equivalent of this data type
tp=type_of(data);
if [tp == 0]
	fprintf('Unknown or unexpected data type: %s\n',class(data));
	if strcmp(class(data),'cell')
		if iscellstr(data)
			tp=7;
		else
			data=-1;
			return;
		end
	else
		data=-1;
		return;
	end
end

% handle each class of data differently
switch tp

	% 'normal', non-complex data types
	case {1,2,3,4,5,10,11,12,13,14,15},
		si=size(data);
		sisi=size(si);
		hdr=[uint64(tp),uint64(sisi(2)),uint64(si)];
		dstd=typecast(hdr,'uint8');
		if strcmp(class(data),'logical')
			data=uint8(data);
		end
		dstd=[dstd,typecast(reshape(data,[1,numel(data)]),'uint8')];

	% complex data types
	case {6,9},
		si=size(data);
		sisi=size(si);
		hdr=[uint64(tp),uint64(sisi(2)),uint64(si)];
		dstd=typecast(hdr,'uint8');
		junk=zeros(1,2,n_elements(data),class(data));
		junk(1,1,:)=real(data(:));
		junk(1,2,:)=imag(data(:));
		dstd=[dstd,typecast(reshape(junk,[1,numel(junk)]),'uint8')];

	% string
	case 7,

		% a single string is of type char
		if strcmp(class(data),'char')
			hdr=[uint64(7),uint64(1),uint64(1)];
			dstd=typecast(hdr,'uint8');
			dstd=[dstd,typecast(uint32(numel(data)),'uint8')];
			dstd=[dstd,uint8(data)];

		% string arrays are type cell
		else

			% prepare the header - it contains the number of strings
			si=size(data);
			sisi=size(si);
			hdr=[uint64(tp),uint64(sisi(2)),uint64(si)];
			dstd=typecast(hdr,'uint8');

			% append the length of each element
			for i=1:numel(data);
				str=char(data(i));
				si=size(str);
				dstd=[dstd,typecast(uint32(si(2)),'uint8')];
			end

			% append the actual string
			for i=1:numel(data),
				str=char(data(i));
				dstd=[dstd,uint8(str)];
			end
		end

	% structure
	case 8,

		% the header contains the size of the structure
		si=size(data);
		sisi=size(si);
		hdr=[uint64(tp),uint64(sisi(2)),uint64(si)];
		dstd=zeros(1,(numel(hdr)*8),'uint8');
		dstd(1,1:numel(hdr)*8)=typecast(hdr,'uint8');

		% it is followed by the tags
		tags=fieldnames(data);
		dstd=[dstd,typecast(uint32(numel(tags)),'uint8')];
		cdat=std_data_fmt(tags);
		dstd=[dstd,cdat];

		% add the actual data
		for i=1:numel(data),
			for j=1:numel(tags),
				str=char(tags(j));
				cdat=getfield(data(i),str);
				dstd=[dstd,std_data_fmt(cdat)];
			end
		end

	% all other cases - we should never get here
	otherwise,
		fprintf('Unknown or unexpected data type: %s\n',class(data));
		dstd=-1;
		return;
end
