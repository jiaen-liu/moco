function y=im_intp_res(im,ratio_res,varargin)
    if all(ratio_res==1)
        y=im;
        return;
    end
    si=size(im);
    nx=si(1);
    ny=si(2);
    nv=1;
    if numel(si)==2
        nz=1;
    elseif numel(si)==4
        nz=si(3);
        nv=si(4);
    elseif numel(si)==3
        nz=si(3);
    end
    if numel(ratio_res)==2
        ratio_res(3)=1;
    end

    method='linear';
    nvarargin=length(varargin);
    if nvarargin>0 && ~isempty(varargin{1})
        method=varargin{1};
    end
    
    im=reshape(im,[nx,ny,nz,nv]);
    xo=[0:nx-1]-(nx-1)/2;
    yo=[0:ny-1]-(ny-1)/2;
    zo=[0:nz-1]-(nz-1)/2;

    nxn=round(nx*1/ratio_res(1));
    nyn=round(ny*1/ratio_res(2));
    nzn=round(nz*1/ratio_res(3));

    xn=([0:nxn-1]-(nxn-1)/2)*ratio_res(1);
    yn=([0:nyn-1]-(nyn-1)/2)*ratio_res(2);
    zn=([0:nzn-1]-(nzn-1)/2)*ratio_res(3);

    [Xo,Yo,Zo]=ndgrid(xo,yo,zo);
    [Xn,Yn,Zn]=ndgrid(xn,yn,zn);

    co=cat(4,Xo,Yo,Zo);
    cn=cat(4,Xn,Yn,Zn);

    y=zeros(nxn,nyn,nzn,nv);
    for i=1:nv
        y(:,:,:,i)=interp3_nmat(co,cn,reshape(double(im(:,:,:,i)),[nx,ny,nz]),method);
    end
end