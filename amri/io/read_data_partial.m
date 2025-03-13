function data=read_data_partial(file,doff,partial)

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
%     October 8, 2020
%     Jacco de Zwart
%     LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%     E-mail: Jacco.deZwart@nih.gov
%
% MODIFICATION HISTORY:
%     #1    Small change to the error message when testing on partial(2) fails.
%           JAdZ, October 9, 2020
%-

% partial should be a 3-element array
if [ numel(partial) ~= 3 ]
	fprintf('ERROR: Partial should contain 3 elements!\n');
	return;
end

% open the file
lun=fopen(file,'r');
if [ lun == -1 ]
	fprintf('Error opening file %s\n',file);
	data=-1;
	return;
end

% read the data descriptor
fseek(lun,doff,'bof');
dhdr=fread(lun,2,'int64=>int64');
fseek(lun,-16,'cof');
dhdr=fread(lun,dhdr(2)+2,'int64=>int64');
dtyp=dhdr(1);
ndim=dhdr(2);
dims=dhdr(3:end);

% check validity of partial
if [ ndim < partial(1) ]
	fprintf('ERROR: Dimension %d does not exist for this %d-dimensional array!\n',partial(1),ndim);
	fclose(lun);
	return;
end
if [ partial(2) < 1 ] || [ partial(2) > dims(partial(1)) ]
	fclose(lun);
	error('ERROR: Beginning index out of range 1..%d!\n',dims(partial(1)));
end
if [ partial(3) < partial(2) ] || [ partial(3) > dims(partial(1)) ]
	fclose(lun);
	error('ERROR: Ending index out of range %d..%d!\n',partial(3),dims(partial(1)));
end

% create data array
pdim=dims;
pdim(partial(1))=partial(3)-partial(2)+1;
udim=int64([0 0 0]);
if [ partial(1) == 1 ]
	udim(1)=1;
else
	udim(1)=dims(1);
	for i=2:(partial(1)-1)
		udim(1)=udim(1)*dims(i);
	end
end
udim(2)=pdim(partial(1));
if [ partial(1) == ndim ]
	udim(3)=1;
else
	udim(3)=dims(partial(1)+1);
	for i=(partial(1)+2):ndim
		udim(3)=udim(3)*dims(i);
	end
end
csiz=udim(1)*type_size(dtyp);
if [ dtyp == 6 ] || [ dtyp == 9 ]
	udim(1)=udim(1)*2;
end
data=zeros(udim,type2class(dtyp));
step=csiz*(dims(partial(1))-udim(2));
skip=csiz*(partial(2)-1);

% read data
if [ skip > 0 ]
	fseek(lun,skip,'cof');
end
for i=1:udim(3)
	data(:,:,i)=fread(lun,[udim(1) udim(2)],type2class(dtyp,'fread_raw',1));
	if [ i < udim(3) ]
		fseek(lun,step,'cof');
	end
end

% make complex if needed
if [ dtyp == 6 ] || [ dtyp == 9 ]
	data=complex(data(1:2:end,:,:),data(2:2:end,:,:));
end

% reshape
pdim=reshape(pdim,[1 ndim]);
data=reshape(data,pdim);

% close file
fclose(lun);
