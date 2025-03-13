function ste_info=ste_motion(ste_info,sorted_ste,par)
% interpolate the clustering index for each shot
    idx_shot_cenpf=sorted_ste.idx_shot_cenpf;
    idx_shot_cenff=sorted_ste.idx_shot_cenff;
    idx=ste_info.idx;
    mot_par=ste_info.mot_par;
    [~,nv_cyc,nvf]=size(mot_par);
    if isfield(par,'rect_sliding_window') && par.rect_sliding_window>0
        for i=1:6
            mot_par(i,:,:)=ste_rect_sliding_window(reshape(mot_par(i,:,:),...
                                                           [nv_cyc,nvf]),...
                                                   par.rect_sliding_window);
        end
    end
    para=sorted_ste.para;
    if isfield(par,'t_mot_filt') && par.t_mot_filt>0
        % low pass filter the motion curve
        fs=1/(sorted_ste.nshot_v*para.tr*1e-3);
        for i=1:6
            mot_filted=butter_fil(a2v(mot_par(i,:,:)),...
                                  1/par.t_mot_filt,...
                                  fs,1,6,'low');
            mot_par(i,:,:)=reshape(mot_filted,[nv_cyc,nvf]);
        end
    end
    nv=sorted_ste.nv;
    nshot=sorted_ste.nshot;
    nk_shot=para.nk_shot;
    ste_intp_res=par.ste_intp_res;
    s2=para.sense_rate_s;
    s1=para.sense_rate_p;
    dkz_caipi=para.dkz_caipi;
    n_interl=para.n_interleaves;
    nc=ste_info.nc;
    fcrptd=ste_info.fcrptd;
    idx_gre=ste_info.idx_gre;
    % 
    if ~sorted_ste.combined
        if numel(idx_shot_cenff)<2
            fcrptd_gre=false(nshot,1);
        else
            fcrptd_gre=interp1(idx_shot_cenff,...
                               double(fcrptd(1:numel(idx_shot_cenff))),...
                               [1:nshot],...
                               'nearest', ...
                               'extrap')>0.5;
        end
    else
        fcrptd_gre=fcrptd(sorted_ste.idx_shot(1,:));
    end

    fcrptd_gre=fcrptd_gre(:);
    mot_par_tmp=reshape(mot_par,[6,nv_cyc*nvf]).';
    % interpolate the motion parameters
    if ~sorted_ste.combined && ~field_true(par,'cluster_b0')
        mot_par_gre=interp1(idx_shot_cenpf,mot_par_tmp(1:nv,:),[1:nshot],'pchip').';
    else
        mot_par_gre=mot_par_tmp(sorted_ste.idx_shot(2,:),:).';
    end
    mot_par_gre(1:3,:)=mot_par_gre(1:3,:); % radian
    mot_par_gre(4:6,:)=mot_par_gre(4:6,:)*ste_intp_res*1e-3; % meter
    % nk of shot x n shot
    if isfield(sorted_ste,'kyz') && ~isempty(sorted_ste.kyz)
        ky=combine_dim(sorted_ste.kyz(1,1:para.nk_shot,:),...
                      [1,2]);
        kz=combine_dim(sorted_ste.kyz(2,1:para.nk_shot,:),...
                      [1,2]);
    else
        if para.isgre
            kyz_tmp=gen_pe_sense(para.np,para.n_partitions,...
                                 s1,s2,dkz_caipi,0,0);
            ky=kyz_tmp(1,:);
            kz=kyz_tmp(2,:);
        else
            i_interl=mod([0:nshot-1],n_interl);
            ik_shot=[0:nk_shot-1].';
            nshot_rep=n_interl*para.n_partitions/s2;
            i_part=floor(mod([0:nshot-1],nshot_rep)/n_interl)*s2-floor(para.n_partitions/2);
            ky=i_interl*s1+ik_shot*s1*n_interl-floor(para.np/2);
            kz=mod(dkz_caipi*(i_interl+ik_shot*n_interl),s2)+i_part;
        end    
    end
    cshot=zeros(nc,1);
    % store motion, k coordinate information
    md=cell(nc,1);

    for i=1:nc
        % assign motion to clusters
        idx_gre_cl=find(idx_gre==i);
        md{i}.m=reshape(mot_par_gre(:,idx_gre_cl),...
                        [3,2,numel(idx_gre_cl)]);
        % assign k coordinates to clusters
        ky_tmp=ky(:,idx_gre_cl);
        kz_tmp=kz(:,idx_gre_cl);
        md{i}.k=[ky_tmp(:),kz_tmp(:)].';
        md{i}.nk=numel(idx_gre_cl)*nk_shot;
        % find the noncorrupted shots
        idx_gre_cl_noncrpt=find(idx_gre==i&(~fcrptd_gre));
        if isempty(idx_gre_cl_noncrpt)
            idx_gre_cl_noncrpt=idx_gre_cl;
        end
        [~,ic]=min(squeeze(min((ky(:,idx_gre_cl_noncrpt).^2+kz(:,idx_gre_cl_noncrpt).^2).^0.5,[],1)));
        cshot(i)=idx_gre_cl_noncrpt(ic);
        % find the full-fov whose position is closest to the shot
        % to get sensitivity
        [~,ii]=min(abs(idx_shot_cenpf-idx_gre_cl(ic)));
        % md{i}.isen=find_nearest_fullfov(ii,nv_cyc,fcrptd);
        if field_true(par,'test_sense')
            % make sure the full nav for the center is not corrupted.
            if fcrptd(floor((ii-1)/nv_cyc)+1)
                warning('*** This is in test_sense mode. The volume was corrupted! ***');
            end
            md{i}.isen=floor((ii-1)/nv_cyc)+1;
        else
            if strcmp(par.vendor,'philips') && ...
                    ~((sorted_ste.y_cyc~=0) && (sorted_ste.z_cyc~=0))
                md{i}.isen=1;
            else
                md{i}.isen=find_closest_pos(mot_par,mod(ii-1,nv_cyc)+1,floor((ii-1)/nv_cyc)+1,fcrptd);
            end
        end
        
        
    end
    ste_info=addvar2struct(ste_info,...
                           'cshot','md');
end
