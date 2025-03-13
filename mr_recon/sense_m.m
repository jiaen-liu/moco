% 20240319: introduce local polynomial fitting to the sensitivity data
function [sensit,ref,cp,imchan,cp_im]=sense_m(data,mask,cov_mat,varargin)

    p=inputParser;
    % apply same phase combination to all volumes
    % otherwise, each volume has its own phase combination
    cp_all=0;
    ref_in=[];
    ref_pha_in=[];
    ress=4;
    frac_sense_sm=8;
    cp_z_shift=0;
    poly_fit=0;
    % unit: mm, a quater wave length
    % 7T: 30 mm
    % 10.5T: 20 mm
    poly_fit_ker=20; 
    addParameter(p,'cp_all',cp_all,@isnumeric);
    addParameter(p,'ref_in',ref_in,@isnumeric);
    addParameter(p,'ref_pha_in',ref_pha_in,@isnumeric);
    addParameter(p,'ress',ress,@isnumeric);
    addParameter(p,'frac_sense_sm',frac_sense_sm,@isnumeric);
    addParameter(p,'cp_z_shift',cp_z_shift,@isnumeric);
    addParameter(p,'poly_fit',poly_fit,@isnumeric);
    addParameter(p,'poly_fit_ker',poly_fit_ker,@isnumeric);
    p.parse(varargin{:});
    cp_all=p.Results.cp_all;
    ref_in=p.Results.ref_in;
    ref_pha_in=p.Results.ref_pha_in;
    ress=p.Results.ress;
    frac_sense_sm=p.Results.frac_sense_sm;
    cp_z_shift=p.Results.cp_z_shift;
    poly_fit=p.Results.poly_fit;
    poly_fit_ker=p.Results.poly_fit_ker;
    poly_fit_ker=round_odd(poly_fit_ker/ress);
    % 1.414 to be consistent with Peter
    frac=ones(1,3)*frac_sense_sm*1.414;
    mbt=0.2;
    
    [nx,ny,nz,nch,nv]=size(data);
    
    if nz>1
        zext=floor(2*frac(3));
    else
        zext=0;
    end
    
    
    cp=zeros(nch,nv);

    sensit=zeros(nx,ny,nz,nch,nv);
    ref=zeros(nx,ny,nz,nv);

    % normalize data
    data=covNorm(data,cov_mat,4);

    psize=floor(1+0.024*nx);
    % psize=floor(1+0.1*nx);
    cp_im=[];    
    for iv=1:nv
        if size(mask,4)==nv
            mask_iv=(mask(:,:,:,iv));
        else
            mask_iv=mask;
        end
        if total(mask_iv)==0
            continue;
        end

        if isempty(ref_in)
            % calculate reference image magnitude
            ref(:,:,:,iv)=ref_mag(data(:,:,:,:,iv),mask,frac_sense_sm);
            if isempty(ref_pha_in)
                % calculate reference phase
                if cp_all==0 || (cp_all==1 && iv==1)
                    % find the center of the mask
                    nm=total(mask_iv);
                    [xx,yy,zz]=ndgrid([1:nx],[1:ny],[1:nz]);
                    ixc=round(total(xx.*mask_iv)/nm);
                    iyc=round(total(yy.*mask_iv)/nm);
                    mi=zeros(nz,1);
                    center_ph=zeros(nch,nz);
                    % only use central slices
                    % to avoid dark hole problem in 
                    % sensitivity maps
                    % assumping z resolution is 4 mm by default
                    % the selected slab should not be more than 40 mm
                    nz_cp=max(min(ceil(40/ress),nz),1);
                    idx_zcp=idx_truncate(nz,nz_cp)+cp_z_shift;
                    % for iz=floor(nz/4)+1:floor(nz*3/4)
                    % for iz=floor(nz*3/8)+1:floor(nz*5/8)
                    cp_im=zeros(psize*2+1,psize*2+1,nz_cp,nch);
                    cp_im=data(ixc-psize:ixc+psize,iyc-psize:iyc+psize,idx_zcp,:);
                    for iz=idx_zcp
                        mi(iz)=total(mask_iv(ixc-psize:ixc+psize,...
                                             iyc-psize:iyc+psize,...
                                             iz))>0;
                        if mi(iz)
                            imr0=data(ixc-psize:ixc+psize,...
                                      iyc-psize:iyc+psize,...
                                      iz,1);
                            imr0(find(imr0==0))=1;
                            for ich=1:nch
                                center_ph(ich,iz)=angle(total(data(ixc-psize:ixc+psize,...
                                                                   iyc-psize:iyc+psize,...
                                                                   iz,ich)./imr0));
                            end
                        end
                    end
                    for ich=1:nch
                        cp(ich,iv)=angle_ave(center_ph(ich,find(mi)));
                    end
                else
                    cp(:,iv)=cp(:,1);
                end
                ref_c=zeros(nx,ny,nz);
                for ich=1:nch
                    ref_c=ref_c+data(:,:,:,ich,iv).*abs(data(:,:,:,ich,iv))*exp(-1i*cp(ich,iv));
                end
            else
                ref_c=ref_pha_in(:,:,:,iv);
            end
            % combine magnitude and phase
            ref(:,:,:,iv)=ref(:,:,:,iv).*exp(1i*angle(ref_c));
        else
            ref(:,:,:,iv)=ref_in(:,:,:,iv);
        end
        sigx=nx/frac(1);
        sigy=ny/frac(2);
        if nz>1
            win=ndgauss([nx,ny,nz+zext*2],[sigx,sigy,(nz+zext*2)/frac(3)]);
            win=win-win(1,1,1);
            win=win/abs(max(win(:)));
            mask_ext=zeros(nx,ny,nz+2*zext);
            mask_ext(:,:,zext+1:zext+nz)=mask_iv;
            mask_ext(:,:,1:zext)=repmat(mask_iv(:,:,1),[1,1,zext]);
            mask_ext(:,:,1+zext+nz:end)=repmat(mask_iv(:,:,end),[1,1,zext]);
            n=nx*ny*(nz+2*zext);
            mb=abs(fftmr(fftmr(mask_ext,-1,[1,2,3]).*win,1,[1,2,3]));
            mb=mb(:,:,zext+1:nz+zext);
        else
            win=ndgauss([nx,ny],[sigx,sigy]);
            win=win-win(1,1);
            win=win/abs(max(win(:)));
            n=nx*ny;
            mb=abs(fftmr(fftmr(mask_iv,-1,[1,2]).*win,1,[1,2]));
        end
        indx=mb>mbt;
        mb_inv=zeros(nx,ny,nz);
        mb_inv(indx)=1./mb(indx);
        
        mask_iv=mask_iv & abs(ref(:,:,:,iv))~=0;
        mask_fill=imfill(mask_iv, 'holes');
        imchan=zeros(nx,ny,nz+2*zext,nch);
        
        for ich=1:nch
            tmp=data(:,:,:,ich,iv).*mask_iv./...
                (mask_iv.*ref(:,:,:,iv)+1-mask_iv);
            
            if poly_fit
                tmp=sense_poly_fit(tmp,mask_iv,...
                                   (abs(ref(:,:,:,iv))/mean(col(abs(ref(:,:,:,iv))))),...
                                   [poly_fit_ker,poly_fit_ker,1],[1,1,1],2);
                tmp=medfilt3(real(tmp))+1i*medfilt3(imag(tmp));
                tmp=tmp.*mask_fill;
            end
            imchan(:,:,zext+1:nz+zext,ich)=tmp;

            if zext>1
                imchan(:,:,1:zext,ich)=repmat(imchan(:,:,1+zext),[1,1,zext]);
                imchan(:,:,nz+zext+1:end,ich)=repmat(imchan(:,:,nz+zext),[1,1,zext]);
            end
            if nz>1
                schan=fftmr(fftmr(imchan(:,:,:,ich),-1,[1,2,3]).*win,1,[1,2,3]);
            else
                schan=fftmr(fftmr(imchan(:,:,:,ich),-1,[1,2]).*win,1,[1,2]);
            end
            sensit(:,:,:,ich,iv)=schan(:,:,zext+1:zext+nz).*mb_inv;
        end
    end
    
end