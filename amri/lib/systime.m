function date=systime(flag)

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
%     #1    Switched to using datestr, which also works in octave (and on non-
%           *nix platforms).
%           JAdZ, December 8, 2017
%-

% if flag is not defined it is assumed to be zero
if [ nargin == 0 ]
	flag=0;
end

% return current date
junk=datestr(now,0);
junk=datenum(junk);
if [ flag == 1 ]
	epoch=datenum(1970,1,1,0,0,0);
	date=(junk-epoch)*24*60*60;
else
	date=sprintf('%s %s %s %s %s',datestr(junk,8),datestr(junk,3),datestr(junk,7),datestr(junk,13),datestr(junk,10));
end
