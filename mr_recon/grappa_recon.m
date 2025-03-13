function im_recon=grappa_recon(d,kyz,w,par)
    getvar_struct(par,...
                  'nd','ry','rz','dz',...
                  'nkx','nky','nkz');
    nkyz=nky*nkz;
    [ker,ind_uniq]=grappa_ker(ry,rz,dz,nky,nkz);
    nr=length(ker);
    [ind_uniq_item,ind_uniq_a,ind_uniq_b]=unique(ind_uniq);
    nr_uniq=numel(ind_uniq_a);
    [d,kyz]=grappa_format_data(d,kyz,par);
    [nxp,nyp,nzp,nch]=size(d);
    d=reshape(d,[nxp,nyp*nzp,nch]);
    np=nxp*nyp*nzp;
    ikx=[0:nkx-1]-floor(nkx/2);
    Ad=zeros(nxp,nyp,nzp,nkx,nkyz,nch,nr);
    nk_ch=nkx*nkyz*nch;
    % assemble data matrix Ad for convolution with weights
    for izp=1:nzp
        for iyp=1:nyp
            for ixp=1:nxp
                for ir=1:nr
                    %ir=ind_uniq_a(ir_u);
                    indyz=iwrapToN(iyp+ker{ir}.ysn,nyp)+...
                          (iwrapToN(izp+ker{ir}.zsn,nzp)-1)*nyp;
                    Ad(ixp,iyp,izp,:,:,:,ir)=...
                        d(iwrapToN(ixp+ikx,nxp),...
                          indyz,:);
                end
            end
        end
    end
    Ad=reshape(Ad,[np,nk_ch,nr]);
    d=reshape(d,[nxp,nyp,nzp,nch]);
    nx=nxp;
    ny=nyp*ry;
    nz=nzp*rz;
    n=nx*ny*nz;
    im_recon_tmp=zeros(nx*nyp*nzp,nch,nr);
    % apply the convolution
    for ir=1:nr
        % ir_uniq=ind_uniq_b(ir)
        im_recon_tmp(:,:,ir)=Ad(:,:,ir)*w(:,:,ir);
    end
    
    im_recon_tmp=reshape(im_recon_tmp,[nx,nyp,nzp,nch,nr]);
    im_recon_tmp=cat(5,d,im_recon_tmp);
    im_recon_tmp=reshape(im_recon_tmp,[nx,nyp,nzp,nch,ry,rz]);
    im_recon_tmp=permute(im_recon_tmp,[1,5,6,4,2,3]);
    % assign to the right position in k-space
    im_recon=zeros(nx,ny,nz,nch);
    ikyf_start=-floor(ny/2);
    ikzf_start=-floor(nz/2);
    
    for izp=1:nzp
        for iyp=1:nyp
            iyf=mod(kyz(1,iyp,izp)-ikyf_start,ny)+1;
            izf=mod(kyz(2,iyp,izp)-ikzf_start,nz)+1;
% $$$             iyzf=iwrapToN(iyf:iyf-1+ry,ny)+ny*...
% $$$                  (iwrapToN(izf:izf-1+rz,nz)-1);
            im_recon(:,iwrapToN(iyf:iyf-1+ry,ny),...
                     iwrapToN(izf:izf-1+rz,nz),:)=...
                im_recon_tmp(:,:,:,:,iyp,izp);
        end
    end
    
    if nd==3
        im_recon=fftmr(im_recon,-1,[1,2,3])*(nx*ny*nz);
    else
        im_recon=fftmr(im_recon,-1,[1,2])*(nx*ny);
    end
end