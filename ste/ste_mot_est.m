function [mot_par,im_coreg]=ste_mot_est(im,iref,mask,method,fast,res_ste)



    if ~isequalfp(res_ste,res_ste(1),1e-4)
        % interpolate the images to be isotropic
        ste_intp_res=min(res_ste(:));
        im=ste_interp(abs(im),...
                      res_ste,...
                      ste_intp_res*ones(1,3));
        mask=ste_interp(single(mask),res_ste,ste_intp_res*ones(1,3));
        mask=mask>0.5;
        if ~isscalar(iref)
            iref=ste_interp(iref,...
                            res_ste,...
                            ste_intp_res*ones(1,3));
        end
    end

    si=size(im);
    if fast
        nvf=si(5);
        nr=si(4);
        if isscalar(iref)
            im_ref=abs(im(:,:,:,:,iref));
        else
            im_ref=abs(iref);
            if size(im_ref,4)==1
                im_ref=repmat(im_ref,[1,1,1,nr]);
            end
        end
    else
        if numel(si)<4
            nvf=1;
        else
            nvf=si(4);
        end
        nr=1;
        if isscalar(iref)
            im_ref=abs(im(:,:,:,iref));
        else
            im_ref=abs(iref);
        end
    end
    nv=nr*nvf;
    tra_matrix=zeros(4,4,nv);
    mot_par=zeros(6,nv);
    % im_mot_reg=reshape(im_mot_reg,[nxste,nyste,nzste,nvfast]);
    im_coreg=zeros([si(1:3),nr,nvf]);
    tmp1=zeros(4,4,nr,nvf);
    tmp2=zeros(6,nr,nvf);
    
    if fast
        m=feature('numcores');
    else
        m=1;
    end

    im=reshape(im,[si(1:3),nr,nvf]);
    if m>1
        parfor (i=1:nr,m)
            [im_coreg(:,:,:,i,:),tmp1(:,:,i,:),tmp2(:,i,:)]=reg3dv(...
                squeeze(im(:,:,:,i,:)),...
                im_ref(:,:,:,i),...
                'method',method,...
                'mask',mask,...
                'refine',1,...
                'fid',i);
        end
    else
        for i=1:nr
            [im_coreg(:,:,:,i,:),tmp1(:,:,i,:),tmp2(:,i,:)]=reg3dv(...
                squeeze(im(:,:,:,i,:)),...
                im_ref(:,:,:,i),...
                'method',method,...
                'mask',mask,...
                'refine',1,...
                'fid',i);
        end
    end
    if numel(si)>4
        im_coreg=combine_dim(im_coreg,[4,5]);
        tmp1=combine_dim(tmp1,[3,4]);
        tmp2=combine_dim(tmp2,[2,3]);
    end
    
    
    if strcmp(method,'amri')
        mot_par=tmp2;
    else
        for iv=1:nv
            mot_par(:,iv)=afmat2motpar(tmp1(:,:,iv),...
                                       method,...
                                       si(1:3));
        end
    end

    if fast
        mot_par=reshape(mot_par,[6,nr,nvf]);
    end
end
