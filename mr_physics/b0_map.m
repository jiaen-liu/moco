function [b0,pha0]=b0_map(x,dte,pha0)
% dte in second
    if nargout>1 && numel(dte)==1
        error('*** If the phase at TE==0 is requested, please provide the actual echo times! ***');
    end
    
    if nargin>2
        pha0_avail=1;
    else
        pha0_avail=0;
    end
    if pha0_avail
        x=x./exp(1i*pha0);
    end
    dte=dte(:);
    si=size(x);
    nx=si(1);
    ny=si(2);
    necho=si(end);
    if numel(si)>=4
        nz=si(3);
    elseif numel(si)==3
        nz=1;
    end
    x=reshape(x,[nx,ny,nz,necho]);
    mag=abs(x);
    p=angle(x);
    pd=angle(x(:,:,:,2)./x(:,:,:,1));
    pdfit=pd;
    if ~(nx==1 && ny==1)
        pdfit=unwrapper_3d_mask(pdfit);
        pdave=mean(col(pdfit(floor(nx/2)-3:floor(nx/2)+3,...
                             floor(ny/2)-3:floor(ny/2)+3,...
                             floor((nz+1)/2))));
        pdfit=pdfit-(pdave-wrapToPi(pdave));
    end
    b0=zeros(nx,ny,nz);
    if numel(dte)==1
        te=[0:necho-1].'*dte;
    else
        te=dte;
    end
    pp=p(:,:,:,1);
    if ~(nx==1 && ny==1)
        pp=unwrapper_3d_mask(pp);
    end
    % pp=p(:,:,:,1);
    if ~pha0_avail
        pha0=zeros(nx,ny,nz);
    end
    te_norm=(te-te(1))/(te(2)-te(1));
    te_null=te-te(1);
    dte=te(2)-te(1);
    n=nx*ny*nz;
    mag=reshape(mag,[n,necho]);
    p=reshape(p,[n,necho]);
    pdfit=pdfit(:);
    pp=pp(:);
    x=reshape(x,[n,necho]);
    b0_vox=0;
    pha0_vox=0;
    for i=1:n
        mag_voxel = a2v(mag(i, :));
        Atest = [te_null.*mag_voxel, mag_voxel];
        rank_defi=0;
        if rank(Atest)<2
            rank_defi=1;
        end
        if ~any(abs(mag_voxel)==0) && ~rank_defi
            if ~pha0_avail
                dp = wrapToPi(a2v(p(i,:))-...
                              pdfit(i)*te_norm-...
                              pp(i));
                b = dp.*mag_voxel;                
                % A = [te_null.*mag_voxel, mag_voxel];
                A = [te_null.*mag_voxel, mag_voxel];
% $$$                 lastwarn('', '');
                k = A\b;
% $$$                 [warnMsg, warnId] = lastwarn();
% $$$                 if contains(warnMsg,'singular')
% $$$                     disp('*** Caught warning singular! ***');
% $$$                 end
                slope = k(1);
                b0_vox=real(pdfit(i)/2/pi/dte+slope/2/pi);
                pha0_vox=real(fit_pha0(b0_vox,te,col(x(i,:))));
                % repeat
                dp = wrapToPi(a2v(p(i,:))-...
                              b0_vox*te*2*pi-pha0_vox);
                b = dp.*mag_voxel;                
                k = A\b;
                slope = k(1);
                b0_vox=real(b0_vox+slope/2/pi);
                pha0_vox=real(fit_pha0(b0_vox,te,col(x(i,:))));
                
                b0(i) = b0_vox;
                pha0(i)=pha0_vox;
            else
                dp=wrapToPi(a2v(p(i,:))-...
                            pdfit(i)/dte*te);
                b=dp.*mag_voxel;
                A = [te.*mag_voxel];
                k = A\b;
                slope = k(1);
                b0(i) = real(pdfit(i)/2/pi/dte+slope/2.0/pi);

            end
        end
    end
    if ~pha0_avail
        pha0=unwrapper_3d_mask(pha0);
    end
end