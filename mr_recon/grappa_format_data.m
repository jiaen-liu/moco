function [do,kyzo]=grappa_format_data(di,kyz,par)

    getvar_struct(par,'nx','ny','nz','nd',...
                      'ry','rz','dz');

    nyp=ny/ry;
    nzp=nz/rz;
    nyzp=nyp*nzp;
    [~,nyzp_di,nch]=size(di);
    [~,nyzp_k]=size(kyz);
    if nyzp_di~=nyzp || nyzp_k~=nyzp
        error('Input dimension nyz does not match');
    end
    do=zeros(nx,nyp,nzp,nch);
% $$$     iyp_c=ceil(nyp/2);
% $$$     izp_c=ceil(nzp/2);
% $$$     iz_c=ceil(nz/2);
    % find the most central k coordinate
% $$$     [~,ic_kyz]=min(sum(abs(kyz).^2,2));
% $$$     do(:,iyp_c,izp_c,:)=di(:,ic_kyz,:);
    iky_start=-floor(ny/2);
    ikz_start=-floor(nz/2);
    kyzo=zeros(2,nyp,nzp);
    for i=1:nyzp
% $$$         if i==ic_kyz
% $$$             continue;
% $$$         end
        iyp=mod(floor((kyz(1,i)-iky_start)/ry),nyp)+1;
        izp=mod(floor((kyz(2,i)-dz*(iyp-1))/rz),nzp)+1;
        do(:,iyp,izp,:)=di(:,i,:);
        kyzo(:,iyp,izp)=kyz(:,i);
    end
end