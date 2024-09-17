function bytes=size_of(data)

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
%     February 17, 2012
%     Jacco de Zwart
%     LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%     E-mail: Jacco.deZwart@nih.gov
%
% MODIFICATION HISTORY:
%     2018/08/28, JAdZ
%         Added support for structures.
%-

% determine the class
type=class(data);

% handle each class of data differently
switch type

	% single and double have to be tested for being complex
	case {'single','double'}
		bytes=numel(data)*type_size(type);
		if ~isreal(data)
			bytes=bytes*2;
		end

	% cell
	case 'cell'
		numel(data)

	% structure - call recursively
	case 'struct'
		bytes=0;
		tags=fieldnames(data);
		for i=1:numel(tags),
			bytes=bytes+size_of(data(1).(tags{i}));
		end
		bytes=numel(data)*bytes;

	% all other cases
	otherwise
		tsiz=type_size(type);
		if (tsiz == 0)
			fprintf('ERROR: Unable to process data of type "%s"\n',type)
			bytes=-1;
		end
		bytes=numel(data)*tsiz;
end
