function save_data(file,data)

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
%     March 8, 2005
%     Jacco de Zwart
%     LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%     E-mail: Jacco.deZwart@nih.gov
%
% MODIFICATION HISTORY:
%     #1    Now uses typecast in stead of cast.
%           JAdZ, December 16, 2005
%     #2    Can now save complex data.
%           JAdZ, April 2, 2008
%     #3    Bugfix: IP address issue on newer openSUSE systems fixed,
%           which return 127.0.0.2 in addition to the 'real' IP address.
%           JAdZ, February 3, 2010
%     #4    The result of type_of() is checked. If it is zero, an unknown
%           data type, an error is printed.
%           JAdZ, January 31, 2012
%     #5    Support for saving structures. The version number is updated to 11.3
%           for consistency with the current IDL version (11.2 and 11.3 files
%           should be identical).
%           JAdZ, January 31 - February 1, 2012
%     #6    Bugfix: User ID conversion failed in R2018b.
%           JAdZ, September 24, 2018
%     #7    Replaced string() by char() to allow this to work in Octave 5.2.0.
%           JAdZ, March 3, 2020
%     #8    Hostname replaced by string "127.0.0.1", this is obsolete anyway.
%           JAdZ, December 21, 2020
%     #9    No longer storing user ID in the header, it is now always set to
%           30582, aka amri.
%           JAdZ, July 13, 2022
%-

% initialisation
versionnr=11.3;
little_endian=bitget(uint32(1),1);

% create .svd header
head=single(zeros(1,8));
head(1)=single(0);
head(2)=(single(little_endian));
head(3)=single(8);
head(4)=single(versionnr);
%[stat,junk]=unix('id -u');
%junk=char(junk);
%val=str2num(junk)
%head(5)=single(val);
head(5)=30582;
ipaddr=uint8(zeros(1,4));
%[stat,junk]=unix('hostname -i');
%junk=strtok(junk,' ');
junk="127.0.0.1";
for i=1:4,
	[word,junk]=strtok(junk,'.');
	ipaddr(i)=uint8(str2num(word));
end
head(6)=typecast(ipaddr,'single');
head(7:8)=typecast(systime(1),'single');

% call std_data_fmt() to get the 
dbyt=std_data_fmt(data);

% create head2
head2=int64(zeros(1,3));
head2(1)=int64(3);
head2(2)=int64(numel(dbyt));
head2(3)=int64(0);

% open file
fid=fopen(file,'w');
if [ fid == -1 ]
	fprintf('Error opening file %s\n',file);
	data=-1;
	return;
end

% write header, data header, data
fwrite(fid,head,class(head));
fwrite(fid,head2,class(head2));
fwrite(fid,dbyt,'uint8');

% close file
fclose(fid);
