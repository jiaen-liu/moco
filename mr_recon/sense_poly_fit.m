function sfit=sense_poly_fit(s,mask,w,ker_size,res,ord)
% sigma: not used
% larger sigma means increased weighting of distant voxels
% 0.5 to 2
    [nx,ny,nz,nch]=size(s);
    if numel(ker_size)==1
        ker_size=ones(3,1)*ker_size;
    end
    nh=floor(ker_size/2);
    nxl=nx+2*nh(1);
    nyl=ny+2*nh(2);
    nzl=nz+2*nh(3);
    sl=zeros(nxl,nyl,nzl,nch);
    sl(idx_truncate(nxl,nx),...
       idx_truncate(nyl,ny),...
       idx_truncate(nzl,nz),:)=s;
    maskl=false(nxl,nyl,nzl);
    maskl(idx_truncate(nxl,nx),...
          idx_truncate(nyl,ny),...
          idx_truncate(nzl,nz))=mask;
    wl=zeros(nxl,nyl,nzl);
    wl(idx_truncate(nxl,nx),...
          idx_truncate(nyl,ny),...
          idx_truncate(nzl,nz))=w;
    sfit=zeros(nx,ny,nz,nch);
    A=[];
    rx=(-nh(1):nh(1))*res(1);
    ry=(-nh(2):nh(2))*res(2);
    rz=(-nh(3):nh(3))*res(3);
    [rx,ry,rz]=ndgrid(rx,ry,rz);
    
    % wr=exp((-rx.^2/nh(1)^2-ry.^2/nh(2)^2-rz.^2/max(nh(3),1)^2)/sigma^2);
    % generate matrix A to save time
    if ker_size(3)==1
        A=gen_poly(rx(:),ry(:),ord);
    else
        A=gen_poly(rx(:),ry(:),rz(:),ord);
    end
    parfor iz=1:nz
        sl_slice=sl(:,:,iz:iz+2*nh(3),:);
        m_slice=maskl(:,:,iz:iz+2*nh(3));
        w_slice=wl(:,:,iz:iz+2*nh(3));
        sfit_slice=zeros(nx,ny,1,nch);
        for ix=1:nx
            for iy=1:ny
                d=sl_slice(ix:ix+2*nh(1),...
                           iy:iy+2*nh(2),...
                           :,:);
                m=m_slice(ix:ix+2*nh(1),...
                          iy:iy+2*nh(2),:);
                wtmp=w_slice(ix:ix+2*nh(1),...
                             iy:iy+2*nh(2),:);
                % w=abs(d).^2;
                % wtmp=wtmp.*wr;
                [sfit_slice(ix,iy,1,:),~]=sense_poly_fit_ker(d,wtmp,m,ord,res,A);
            end
        end
        sfit(:,:,iz,:)=sfit_slice;
    end
    
end
function [y,A]=sense_poly_fit_ker(d,w,m,ord,res,A)
    [nx,ny,nz,nch]=size(d);
    n=nx*ny*nz;
    d=reshape(d,[n,nch]);
    if numel(w)==n
        w=w(:);
    elseif numel(w)==(n*nch)
        w=reshape(w,[n,nch]);
    end
    if nargin<6
        A=[];
    end
    if isempty(A)
        x=((0:nx-1)-floor(nx/2))*res(1);
        y=((0:ny-1)-floor(ny/2))*res(2);
        z=((0:nz-1)-floor(nz/2))*res(3);
        [x,y,z]=ndgrid(x,y,z);
        x=x(:);
        y=y(:);
        z=z(:);
        if nz==1
            A=gen_poly(x,y,ord);
        else
            A=gen_poly(x,y,z,ord);
        end
        
    end
    % icent=nx*ny*floor(nz/2)+nx*floor(ny/2)+floor(nx/2)+1;
    y=zeros(nch,1);
    if total(m)>=size(A,2)
        for ich=1:nch
            Aw=A(m(:),:).*w(m(:),ich);
            if rank(Aw)<size(Aw,2)
                continue;
            end
            dw=d(m(:),ich).*w(m(:),ich);
            yv=Aw\dw;
            y(ich)=yv(1);
        end
    end
end