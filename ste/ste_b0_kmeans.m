function [idx,c,dist,dist_max,dist_tot,b0_resid,c_b0_cent,c1,c_test,b0_fit,mask,b0_resid_pca]=ste_b0_kmeans(b0,mask,k,ord,res,...
                                                      varargin)
    % k is 
    idx_in=[];
    dist=[];
    dist_max=[];
    dist_tot=[];
    p=inputParser;
    addParameter(p,'idx_in',idx_in,@isnumeric);
    p.parse(varargin{:});
    idx_in=p.Results.idx_in;

    if ~isempty(idx_in)
        k=length(unique(idx_in));
    end
    
    [nx,ny,nz]=size(mask);
    nmsk=total(mask);
    n=nx*ny*nz;
    nv=numel(b0)/n;
    % filter the b0 data
    sfil=5;
    h=fspecial('gauss',sfil,0.6);
    b0=imfilter(b0,h);
    % shrink mask
    hsfil=floor(sfil/2);
    disk=strel('disk',hsfil);
    for i=1:nz
        mask(:,:,i)=imerode(mask(:,:,i),disk);
    end
    % if ord >=6, use the original b0
    if ord>=6
        b0_ori=b0;
    else
        b0_ori=eval_sphere_fitb0(b0,mask,ord,[],res);
    end
    % remove linear component in b0
    [~,b0_resid]=eval_sphere_fitb0(b0,mask,1,[],res);
    [b0,~,~,c_b0]=eval_sphere_fitb0(b0_resid,mask,6,[],res);
    % b0=b0_resid;
    b0_resid=reshape(b0_resid,[n,nv]);
    %
    b0=reshape(b0,[n,nv]);
    b0=b0(mask(:),:);
    b0_resid=b0_resid(mask(:),:);
    nk=length(k);
    dist_tot=zeros(nk,1);
    %% ------- test ------- %%
% $$$     x=([1:nx]-(1+nx)/2)*res(1);
% $$$     y=([1:ny]-(1+ny)/2)*res(2);
% $$$     z=([1:nz]-(1+nz)/2)*res(3);
% $$$     [x,y,z]=ndgrid(x,y,z);
% $$$     A=gen_spher_harm_poly(x(mask),y(mask),z(mask),1);
% $$$     c1n=zeros(4,k(end));
% $$$     % c1n([3],:)=(rand(1,k(end))+2)*21;
% $$$     c1n([2],:)=30;
% $$$     centn=A*c1n;
    %% ------- ---- ------- %%
    for j=1:nk
        if ~isempty(idx_in)
            idx=idx_in;
        else
            [idx,c]=kmeans(b0.',k(j));
            c=c.';
            % calculate distance to centroid for each cluster
            % and the total distance of all samples to their centroid
            dist=cell(k(j),1);
            dist_max=zeros(k(j),1);
            c_b0_cent=zeros((ord+1)^2,k(j));

            for i=1:k(j)
                dist{i}=a2v(mean(abs(b0_resid(:,idx==i)-c(:,i)).^2,1).^0.5);
                dist_max(i)=max(dist{i});
                dist_tot(j)=dist_tot(j)+total(dist{i}.^2);
                %% --- temporal -- %%
                % c_b0_cent will be assigned later
                % c_b0_cent(:,i)=mean(c_b0(:,idx==i),2);
                
            end
            %% ------- test ------- %%
% $$$         c=c+centn;
% $$$         c_b0_cent(1:4,:)=c_b0_cent(1:4,:)+c1n;
            %% ------- ---- ------- %%
            cn=zeros(n,k(j));
            cn(mask(:),:)=c;
            c=reshape(cn,[nx,ny,nz,k(j)]);
        end
        %% --- temporal -- %%
        % there is an issue related to big linear b0
        if j==nk
            c=zeros(nx,ny,nz,k(j));

            for i=1:k(j)
                b0c=mean(b0_ori(:,:,:,idx==i),4);
                [c(:,:,:,i),~,~,c_b0_cent(:,i)]=eval_sphere_fitb0(b0c,mask,ord,[],res);
            end
        end
        %% --- -------- -- %%
    end
    dist_tot=(dist_tot/nv).^0.5;
    b0_resid_tmp=b0_resid;
    b0_resid=zeros(n,nv);
    b0_resid(mask(:),:)=b0_resid_tmp;
    b0_resid=reshape(b0_resid,[nx,ny,nz,nv]);

    % fit each b0 map to the first order relative to the centroid
    c1=zeros(4,nv);
    for i=1:k(end)
        [~,b0_resid_pca,~,c1tmp]=...
            eval_sphere_fitb0(b0_ori(:,:,:,idx==i)-c(:,:,:,i),mask,1,[],res);
        c1(:,idx==i)=c1tmp;
    end
    % adjust the sphere harmonics in the order of x,y,and z
    c1tmp=c1;
    c1(1,:)=c1tmp(1,:)*(0.5*(1/pi)^0.5);
    c1(2,:)=c1tmp(4,:)*(3/4/pi)^0.5;
    c1(3:4,:)=c1tmp(2:3,:)*(3/4/pi)^0.5;

% $$$     c_test=[];
% $$$     b0_fit=[];
% $$$     % validate c_b0_cent is correct
    c_test=zeros(nx,ny,nz,k(end));
    x=([1:nx]-(1+nx)/2)*res(1);
    y=([1:ny]-(1+ny)/2)*res(2);
    z=([1:nz]-(1+nz)/2)*res(3);
    [x,y,z]=ndgrid(x,y,z);
    c_test=reshape(sphere_harm_calc_3d(x,y,z,c_b0_cent),[nx,ny,nz,k(end)]);
    % validate 1st-order b0 fit is correct
    b0_fit=zeros(n,nv);
    for i=1:k(end)
        b0_fit(:,idx==i)=c1(1,idx==i)+c1(2,idx==i).*x(:)+...
            c1(3,idx==i).*y(:)+...
            c1(4,idx==i).*z(:)+...
            a2v(c_test(:,:,:,i));
    end
    b0_fit=reshape(b0_fit,[nx,ny,nz,nv]);
end