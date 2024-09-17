function [data,header]=read_data(file,varargin)

%+
% NAME:
%     
%#
% FUNCTION:
%     Read IDL .svd file
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
%     July 13-20, 2004
%     Jacco de Zwart
%     LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%     E-mail: Jacco.deZwart@nih.gov
%
% MODIFICATION HISTORY:
%     #1    Changed printing to screen from disp(sprintf()) to fprintf().
%           JAdZ, February 28, 2005
%     #2    Changes for largely more efficient memory usage.
%           JAdZ, December 22, 2005
%     #3    Bugfix: Round-off errors for loop counter could cause crash
%           for very small files.
%           JAdZ, January 16, 2007
%     #4    Support for version 11.3 data (no real changes).
%           JAdZ, April 22, 2011
%     #5    Bugfix: One string was printed even when keyword silent was set.
%           JAdZ, February 1, 2012
%     #6    Bugfix: When writing read_mrcamera.m I found that fread should be
%           called with 'type=>type' as the 3rd argument, instead of just
%           'type'. The latter reads data as 'type' but returns the data as
%           double...
%           JAdZ, February 17, 2012
%     #7    All data are now read at once instead of in chuncks. This was
%			originally done because I'd made the same error as under #6 when
%			reading the 'normal' data from disk.
%			JAdZ, April 16, 2012
%     #8    The silent and info parameters are passed on to std_data_rev().
%           JAdZ, August 25, 2015
%     #9    The header is now returned when this is called with [data,header]=.
%           JAdZ, May 25, 2017
%     #10   Keyword partial added to allow reading only a subset of the data.
%           Thanks to Jiaen Liu for the inputParser example code. Optional
%           arguments silent and info are now also processed this way.
%           JAdZ, October 8, 2020
%-

% initialization
max_supported_version=11.3;
little_endian=bitget(uint32(1),1);

% optional arguments 
silent=1;
info=0;
partial=-1;
p=inputParser;
var=[];
%addRequired(p,'file');							% defined above
addOptional(p,'osilent',-1);					% for backward compatibility, where silent was an optional 2nd argument
addOptional(p,'oinfo',-1);						% for backward compatibility, where info was an optional 3rd argument
addParameter(p,'silent',silent,@isnumeric);
addParameter(p,'info',info,@isnumeric);
addParameter(p,'partial',partial,@isnumeric);
p.parse(varargin{:});
if [ p.Results.osilent ~= -1 ] 
	silent=p.Results.osilent;
else
	silent=p.Results.silent;
end
if [ p.Results.oinfo ~= -1 ] 
	info=p.Results.oinfo;
else
	info=p.Results.info;
end
partial=p.Results.partial;

% open file
fid=fopen(file,'r');
if [ fid == -1 ]
	fprintf('Error opening file %s\n',file);
	data=-1;
	return;
end

% first header: the first word indicates if this is a pre-version 9.0 file
junk=fread(fid,1,'float=>float');
if [ junk ~= 0.0 ]
	fprintf('This file was generated with SAVE_DATA version < 9.0.\n');
	fprintf('Files this old are currently not supported.\n');
	return;
end

% first header, 2nd word: little endian y/n = 1/0
junk=fread(fid,1,'float=>float');
if [ junk ~= little_endian ]
	fclose(fid);
	if [ little_endian ]
		fid=fopen(file,'r','ieee-be');
	else
		fid=fopen(file,'r','ieee-le');
	end
	fseek(fid,8,'bof');
end

% first header, 3rd word: length of first header
junk=fread(fid,1,'float=>float');

% read the complete first header
fseek(fid,0,'bof');
head=fread(fid,junk,'float=>float');
if [ head(4) >= 10.0 ] l64=1; else l64=0; end

% second header: long or long64, typically 3 words/12 or 24 bytes, length in 1st element)
if [ l64 ]
	junk=fread(fid,1,'int64=>int64');
	fseek(fid,-8,'cof');
	head2=fread(fid,junk,'int64=>int64');
else
	junk=fread(fid,1,'int32=>int32');
	fseek(fid,-4,'cof');
	head2=fread(fid,junk,'int32=>int32');
end

% check for data newer than supported (use tolerance for testing floating point...)
if [ (head(4)-max_supported_version) > 0.0001 ]
	fprintf('Version conflict: This file contains version %4.1f data, ',head(4));
	fprintf('but this routine supports up to version %4.1f.\n',max_supported_version);
	data=-1;
	return;
end

% show info
if [ info == 1 ]
	fprintf('File info:\n');
	fprintf('  Version: %2.1f\n',head(4));
	fprintf('  64-bit support: %d\n',l64);
	fprintf('  Number of data bytes: %d\n',head2(2));
	fprintf('  Number of header bytes: %d\n',head2(3));
end;

% read data
if [ head2(2) > 0 ]
	if [ numel(partial) == 3 ]
		dstart=ftell(fid);
		dhdr=fread(fid,2,'int64=>int64');
		fseek(fid,-16,'cof');
		dhdr=fread(fid,dhdr(2)+2,'int64=>int64');
		if [ dhdr(1) == 7 ] || [ dhdr(1) == 8 ]
			error('ERROR: Keyword partial not yet implemented for this data type!');	
		end
		if [ dhdr(2) < partial(1) ]
			error('ERROR: Dimension for keyword partial does not exist!');	
		end
		if [ silent ~= 1 ]
			fprintf('Reading file partially: ');
		end
		data=read_data_partial(file,dstart,partial);
		if [ silent ~= 1 ]
			fprintf('done\n');
		end
	else
		if [ partial(1) ~= -1 ]
			fprintf('WARNING: Malformatted keyword partial, reading all of the data!\n');
		end
		if [ silent ~= 1 ]
			fprintf('Reading contents of file:       ');
		end
		data=fread(fid,head2(2),'uint8=>uint8');
		if [ silent ~= 1 ]
			fprintf('\n');
			fprintf('Converting to the proper data format: ');
		end
		data=std_data_rev(data,silent,info);
		if [ silent ~= 1 ]
			fprintf('done\n');
		end
	end
else
	data=-1;
	fclose(fid);
	return;
end

% read header
if [ head2(3) > 0 ]
	header=fread(fid,head2(3),'uint8=>uint8');
	header=std_data_rev(header);
else
	header=-1;
end

% close file
fclose(fid);
