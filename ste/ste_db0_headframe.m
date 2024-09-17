function [db0,im_warp,mask_in_fov,im_recov]=ste_db0_headframe(im,iref,mot,dte,ncontr,res)
    breal=isreal(im);
    mot=squeeze(mot);
    [~,nv]=size(mot);
    si=size(im);
    nx=si(1);
    ny=si(2);
    nz=si(3);
    im=reshape(im,[nx,ny,nz,ncontr,nv]);
    nvim=nv;
    ncontrim=ncontr;
    im_warp=zeros(nx,ny,nz,ncontr,nv);
    if nargout>3
        im_recov=zeros(nx,ny,nz,ncontr,nv);
    end
    % align to head frame
    for i=1:nv
        im_warp(:,:,:,:,i)=move_matrix(im(:,:,:,:,i),...
                                       mot(:,i),...
                                       'inverse','res',res);
        if nargout>3
            im_recov(:,:,:,:,i)=move_matrix(im_warp(:,:,:,:,i),...
                                            mot(:,i),...
                                            'forward','res',res);
        end
    end
    % calculate b0 difference
    db0=zeros(nx,ny,nz,nv);
    if numel(size(iref))>=3
        im_ref=iref;
    else
        im_warp(:,:,:,:,iref)=im(:,:,:,:,iref);
        im_ref=im_warp(:,:,:,:,iref);
    end
    for i=1:nv
        if ncontr==1
            db0(:,:,:,i)=unwrapper_3d_mask(...
                squeeze(angle(im_warp(:,:,:,1,i)./...
                              im_ref(:,:,:,1))),...
                ones(si(1:3)))/2/pi/dte;
        else
            db0(:,:,:,i)=unwrapper_3d_mask(...
                squeeze(angle(im_warp(:,:,:,2,i)./...
                              im_warp(:,:,:,1,i)./...
                              im_ref(:,:,:,2).*...
                              im_ref(:,:,:,1))),...
                ones(si(1:3)))/2/pi/dte;
        end
    end
    mask_in_fov=abs(im_warp(:,:,:,1,2:end-1))>0 & ...
        ~isnan(abs(im_warp(:,:,:,1,2:end-1)));
    mask_in_fov=all(mask_in_fov,5);
end