function [data,bytes]=std_data_rev(data,varargin)

%+
% NAME:
%     
%#
% FUNCTION:
%     Convert standardized data back to an array of the proper type.
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
%     July 15-21, 2004
%     Jacco de Zwart
%     LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%     E-mail: Jacco.deZwart@nih.gov
%
% MODIFICATION HISTORY:
%     #1    Small changes so that it works properly in Matlab 7.0
%           (release 14).
%           JAdZ, July 27, 2004
%     #2    Eliminated for-loop, improving performance about 100-fold.
%           JAdZ, February 28, 2005
%     #3    Now uses typecast in stead of cast.
%           JAdZ, December 16, 2005
%     #4    Small improvements in memory usage efficiency. Fixes to allow
%           reading of Peter's data with the wrong header size.
%           JAdZ, December 23, 2005
%     #5    Can now handle complex and dcomplex data.
%           JAdZ, September 27, 2006
%     #6    Support for files containing structures. The number of bytes in data
%           that was used is now also returned.
%           JAdZ, November 25 & 28, 2011
%     #7    Bugfix: Only the first characer of each structure tag name was
%           retrieved.
%           JAdZ, January 31, 2012
%     #8    Bugfix: Reading a structure containing complex data failed since the
%           size of the complex data element was not computed correctly. Multi-
%           dimensional string arrays are now supported.
%           JAdZ, February 1, 2012
%     #9    Bugfix: There was an erroneous extra first dimension of size 1 when
%           processing string-type data.
%           JAdZ, July 9, 2014
%     #10   Added silent and info optional arguments. Some values converted to
%           unit64 data type.
%           JAdZ, August 25, 2015
%     #11   Fixed an integer type conflict (cannot add int32 to uint64?). String
%           used for structure field names (probably in all cases) needed a
%           transpose to be correct.
%           JAdZ & Catie Chang, August 18, 2016
%     #12   Bugfix: One file offset needed to be explicity converted to a uint64
%           otherwise a crash would occur.
%           JAdZ, May 25, 2017
%     #13   Significant performance increase when evaluating a string version of
%           the command that copies a large part of the data (~line 250).
%           Jiaen Liu, May 21, 2018
%-

% Resolve array size overflow issue, by Jiaen Liu, Aug-19, 2024

% optional arguments
silent=1;
info=0;
if nargin > 1
	silent=varargin{1};
	info=varargin{2};
end

% read header to determine array type and size
idltype=typecast(uint8(reshape(data(1:8),1,8)),'int64');
if idltype > int64(72057594037927936)
	idltype=swap_endian(typecast(uint8(reshape(data(1:8),1,8)),'int64'));
	swap=1;
    if info == 1
        fprintf('Swapping required\n');
    end
else
	swap=0;
end
if idltype > int64(2147483647)
	idltype=typecast(uint8(reshape(data(1:4),1,4)),'int32');
	if swap idltype=swap_endian(idltype); end
	if idltype > 256
		idltype=idltype-256;
		dims=typecast(uint8(reshape(data(5:12),1,8)),'int64');
		if swap dims=swap_endian(dims); end
		dims=double(dims);
		idlsize=int64(zeros(1,dims+3));
		for i=1:dims,
			junk=typecast(uint8(reshape(data(i*8+5:(i+1)*8+4),1,8)),'int64');
			if [ swap ] junk=swap_endian(junk); end
			idlsize(i+1)=junk;
		end
		hdrsize=(dims+1)*8+4;
	else
		dims=typecast(uint8(reshape(data(5:8),1,4)),'int32');
		if [ swap ] dims=swap_endian(dims); end
		dims=double(dims);
		idlsize=int64(zeros(1,dims+3));
		for i=1:dims,
			junk=typecast(uint8(reshape(data((i+1)*4+1:(i+2)*4),1,4)),'int32');
			if [ swap ] junk=swap_endian(junk); end
			idlsize(i+1)=int64(junk);
		end
		hdrsize=(dims+2)*4;
	end
else
	dims=typecast(uint8(reshape(data(9:16),1,8)),'int64');
	if [ swap ] dims=swap_endian(dims); end
	dims=double(dims);
	idlsize=int64(zeros(1,dims+3));
	for i=1:dims,
		junk=typecast(uint8(reshape(data((i+1)*8+1:(i+2)*8),1,8)),'int64');
		if [ swap ] junk=swap_endian(junk); end
		idlsize(i+1)=junk;
	end
	hdrsize=(dims+2)*8;
end
idlsize(1)=int64(dims);
idlsize(dims+2)=idltype;
idlsize(dims+3)=1;
for i=1:dims,
	idlsize(dims+3)=double(idlsize(dims+3))*double(idlsize(i+1));
end
if [ info == 1 ]
    fprintf('\n  Number of dimension:s %d\n',idlsize(1));
       fprintf('    [ ');
    for i=1:dims,
	    fprintf('%d ',idlsize(i+1));
    end
    fprintf(']\n');
    fprintf('  IDL data type: %d\n',idlsize(dims+2));
end

% determine the data type
type=double(idlsize(dims+2));
cls=type2class(type);

% string-type data
if [ type == 7 ]

	% extract the strings
	nstr=idlsize(idlsize(1)+3);
	slen=typecast(data(hdrsize+1:hdrsize+(4*nstr)),'int32');
	offs=int32(hdrsize+(4*nstr));
	for i=1:nstr,
		str=native2unicode(data(offs+1:offs+slen(i)));
		str=str';
		offs=offs+slen(i);
		if (i == 1)
			dat=str;
		else
			dat=char(dat,str);
		end
	end
	dat=cellstr(dat);

	% reshape data if needed
	if [idlsize(1) > 1]
		dat=reshape(dat,idlsize(2:(idlsize(1)+1)));
	end

	% return data and bytes used
	data=dat;
	bytes=offs;

	return;

% struct-type data
elseif [ type == 8 ]

	% determine the number of tags
	ntag=typecast(data(hdrsize+1:hdrsize+4),'int32');

	% extract the tag names and their size in file
	[tags,offs]=std_data_rev(data(hdrsize+5:end));
	offs=uint64(offs);

	% determine offset to the structure data in bytes
	offs=offs+uint64(hdrsize+4);

	% compute the number of structure elements
	nelm=1;
	for i=1:idlsize(1)
		nelm=nelm*idlsize(i+1);
	end

	% read the structure content
	for i=1:nelm,
		for j=1:ntag,
            ii1=offs+1;
            cmd=['[elem,coff]=std_data_rev(data(' int2str(ii1) ':end));'];
            eval(cmd);
			% [elem,coff]=std_data_rev(data(offs+1:end));
			offs=offs+uint64(coff);
			if (i == 1)
				if (j == 1)
					dat=struct(lower(char(tags(j))),elem);
				else
					dat=setfield(dat,lower(char(tags(j))),elem);
				end
			else
				dat(i).(lower(char(tags(j))))=elem;
			end
		end
	end

	% return data and bytes used
	data=dat;
	bytes=offs;
else

	% if these are 'normal' data, convert them
	junk=uint64(1);
	for i=2:(dims+1),
		junk=junk*uint64(idlsize(i));
	end
	junk=junk*uint64(type_size(cls));
	if [ element_of(type,[6,9]) ]
		junk=junk*uint64(2);
	end
	if info == 1
		fprintf('  Data size in bytes (after header stripping): %d\n',junk);
		fprintf('  Word size in bytes: %d\n',type_size(cls));
	end
	% NOTE: THE NEXT LINE CONSUMES TOO MUCH MEMORY! 6x ARRAY SIZE!!!
% $$$     ndata=numel(data);
% $$$     data=data(1:end-(ndata-hdrsize-junk));
% $$$     data=data(hdrsize+1:end);
% $$$ 	ii1=hdrsize+1;
% $$$ 	ii2=hdrsize+junk;
	%data=data(hdrsize+1:hdrsize+junk);
    %% 2024-08-19
    % Resolve array size overflow issue, by Jiaen Liu
    datatmp=zeros(junk,1,'uint8');
    max_seg=8*1024^3;
    n_seg=ceil(double(junk)/max_seg);
    for i_seg=1:n_seg
        ii1_target=(i_seg-1)*max_seg+1;
        ii1_source=ii1_target+hdrsize;
        if i_seg<n_seg
            ii2_target=i_seg*max_seg;
        else
            ii2_target=double(junk);
        end
        ii2_source=ii2_target+hdrsize;
        cmd=['datatmp(',...
             int2str(ii1_target),':',int2str(ii2_target),')',...
             '=data(',...
             int2str(ii1_source),':',int2str(ii2_source),');'];
        eval(cmd);
    end
    data=datatmp;
    clear datatmp;
	if [ element_of(type,[1,2,3,4,5,6,9,12,13,14,15]) ]
		wordsi=uint64(type_size(cls));
		junk=uint64(idlsize(dims+3))*wordsi;
		data=typecast(data,cls);
		if [ element_of(type,[6,9]) ]
			junk=junk*uint64(2);
			data=reshape(data,[2,n_elements(data)/2]);
			data=complex(data(1,:),data(2,:));
		end
		if [ idlsize(1) > 1 ]
			data=reshape(data,double(idlsize(2:(dims+1))));
		end
		if [ swap ] data=swap_endian(data); end

		% compute the number of bytes used
		bytes=hdrsize+junk;
	else
		if [ cls(1) == '' ]
			fprintf('This file contains an unknown data type.\n');
		else
			fprintf('This file contains data of type %s, which is yet not supported in the matlab version.\n',cls);
		end
	end
end
