function im_recon=sense_recon(imd,sensit,mask,inv_cov,s1,s2,ds,dy,dz)
% s1 and s2: sense factor in two dimensions    
% ds: cappi shift
% dy: offset of ky
% dz: offset of kz
% sensit was calculated as
% sensit=sense_m(b1d,mask,eye(nch,nch));
    [nx,ny,nz,nch]=size(sensit);
    n=nx*ny*nz;
    nyz=ny*nz;
    smat=sense_mat(sensit,inv_cov,s1,s2,ds,dy,dz);
    % mask in k space
    mk=sense_mask_k(ny/s1,nz/s2,s1,s2,ds,dy,dz);
    % image data in k space
    [~,nys,nzs,~]=size(imd);
    ns=nx*nys*nzs;
    nyzs=nys*nzs;
    imd=reshape(imd,[ns,nch]);
    imd=imd*conj(inv_cov);
    imd=reshape(imd,[nx,nyzs,nch]);
    b=zeros(nx,nyz,nch);
    b(:,mk(:),:)=imd;
    b=reshape(b,[nx,ny,nz,nch]);
    b=fftmr(b,-1,[1,2,3]);
    b=sum(reshape(conj(sensit).*b,[n,nch]),2);
    im_recon=reshape(smat\b,[nx,ny,nz]);
    im_recon(isnan(im_recon))=0;
    im_recon(~mask)=0;
    mask_blur=zeros(nx,ny,nz);
    h=fspecial('gauss',5,1);
    for i=1:nz
        mask_blur(:,:,i)=imfilter(double(mask(:,:,i)),h);
    end
    im_recon=im_recon.*mask_blur;
end