function [ste_info,b0_fit]=ste_cluster(ste_info,par,...
                                       para,sorted_ste,...
                                       db0,mask_brain)
    b0_fit=[];

    nc0=par.nc0;
    nc_max=par.nc_max;
    max_t_km=par.max_t_km;
    max_a_km=par.max_a_km;
    max_t=par.max_t;
    max_a=par.max_a;
    max_f=par.max_f;
    
    nte=para.nte_contr;
    if isfield(par,'pcaflag')
       pcaflag = par.pcaflag; 
       if pcaflag
           %nc = 1;
           %Number of PCA modes to be used. To be provided, maybe add automatic choice of number of PCA Modes Later
           npca = par.npca;
       end
    else
        pcaflag = 0;
    end
    
    mot_par=ste_info.mot_par;
    [~,nv_cyc,nvf]=size(mot_par);
    nv=sorted_ste.nv;
    fcrptd=ste_info.fcrptd;
    resx_ste=sorted_ste.para.steref_res_r;
    resy_ste=sorted_ste.para.steref_res_p;
    resz_ste=sorted_ste.para.steref_res_s;
    nshot=sorted_ste.nshot;
    idx_shot_cenpf=sorted_ste.idx_shot_cenpf;
    if ~field_true(par,'cluster_b0')
        % cluster the motion parmeters into a few groups
        % within each, there is limited motion and the field
        % changes can be approximated up to 1st order changes.
        % k-means clustering
        % first cluser based on intact full-fov volumns
        % second fill the moved volumns
        if isfield(par,'idx') && ~isempty(par.idx)
            % use predefined parameters in the configuration file

            idx=zeros(nv_cyc*nvf,1);
            idx_par=par.idx(:);
            nc=max(idx_par);
            if length(idx)>length(idx_par)
                idx(1:length(idx_par))=idx_par;
                idx(length(idx_par)+1:end)=idx_par(end);
            elseif length(idx)<length(idx_par)
                idx=idx_par(1:length(idx));
            else
                idx=idx_par;
            end
        else
            nc=nc0-1;
            nvf_noncrpt=total(~fcrptd);
            mot_par_km=reshape(squeeze(mot_par(:,floor(end/2),~fcrptd)),[6,nvf_noncrpt]).';
            for i=1:nc_max-nc0+1
                nc=nc+1;
                if nvf_noncrpt<nc
                    error('There are not enough static full-fov volumns to cluster!');
                end
                [idx_km,c]=kmeans(mot_par_km,nc);
                max_tmp=max(abs(mot_par_km-c(idx_km,:)),[],1);
                max_t_tmp=max(max_tmp(1:3));
                max_a_tmp=max(max_tmp(4:6));
                if max_t_tmp<=max_t_km && max_a_tmp<=max_a_km
                    break;
                end
            end
            idx=zeros(nv_cyc,nvf);
            idx(:,~fcrptd)=repmat(idx_km,[1,nv_cyc]).';
            for i=1:nvf
                if ~fcrptd(i)
                    continue;
                end
                for j=1:nv_cyc
                    idx(j,i)=idx(1,find_closest_pos(mot_par,j,i,fcrptd));
                end
            end
            idx=idx(:);
        end
    else
        % cluster based on b0
        res=[resx_ste,resy_ste,resz_ste];
        db0=combine_dim(db0,[4,5]);
        db0=db0(:,:,:,1:nv);
        if isfield(par,'cluster_b0_auto') && par.cluster_b0_auto>0
            % determine the number of clusters automatically based on b0 changes
            [idx,c,dist,dist_max,dist_tot,~,~,~,~,~]=...
                ste_b0_kmeans(db0,mask_brain,[1:20],par.ord,res);
            dist_tot_bottom=mean(dist_tot(end-2:end));
            nc=find(abs(dist_tot-dist_tot_bottom)<...
                    (par.cluster_b0_auto)*abs(dist_tot(1)-dist_tot_bottom),1,'first');
            disp(['*** ' int2str(nc) ' clusters were automatically determined. ***']);
        else
            nc=par.nc_max;
        end
        if field_true(par,'test_sense')
            idx=ste_idx_ibreak(sorted_ste);
            [~,c,dist,dist_max,dist_tot,~,c_b0_cent,c1,~,b0_fit]=...
                ste_b0_kmeans(db0,mask_brain,[nc],par.ord,res,'idx_in',idx);
        else
            % test
% $$$             idx(1:528)=1;
% $$$             idx(529:1368)=2;
% $$$             idx(1369:1476)=3;
% $$$             idx=idx(:);
            [idx,c,dist,dist_max,dist_tot,~,c_b0_cent,c1,~,b0_fit,...
             mask_brain_pca,b0_resid_pca]=...
                ste_b0_kmeans(db0,mask_brain,[nc],par.ord,res);
        end
        scores=[];
        coeffs=[];
        if pcaflag && nc==1
            [nx,ny,nz,nb0]=size(b0_resid_pca);
            scores=zeros(nx,ny,nz,npca,nte);
            coeffs=zeros(nb0,npca,nte);
            for ie=1:nte
                %define te
                te=para.te_contr(ie);
                %Perform PCA
                [scorestmp,coeffstmp,~,~,b0_fit] = ste_pca_b0(b0_resid_pca,mask_brain_pca,npca,te);
                scores(:,:,:,:,ie) = scorestmp;
                coeffs(:,:,ie) = coeffstmp;
                %other outputs of RunPCAB0 are values that may later be
                %used to pick the number of pca modes automatically.
            end
        elseif pcaflag && nc>1
            error('*** Only one cluster is allowed when PCA is used for B0 correction! ***');
        end
    end
    % generate a report about max within cluster rotation
    mot_cluster=cell(nc,1);
    tmp_mot_par=combine_dim(mot_par,[2,3]);
    for i=1:nc
        idx_cluster=find(idx(:)==i);
        mot_par_cluster=tmp_mot_par(:,idx_cluster);
        mot_par_ave_cluster=mean(mot_par_cluster,2);
        mot_par_max_cluster=max(abs(mot_par_cluster-mot_par_ave_cluster),[],2);
        mot_cluster{i}=var2struct('idx_cluster','mot_par_cluster',...
                                  'mot_par_ave_cluster','mot_par_max_cluster');
    end
% prepare the data struction for gre reconstruction
% total number of shots in the gre data
    idx_gre=zeros(nshot,1);
    % interpolate the clustering index for each shot
    if ~sorted_ste.combined && ~field_true(par,'cluster_b0')
        idx_gre=interp1(idx_shot_cenpf,idx(1:nv),[1:nshot],'nearest', ...
                        'extrap');
        idx_gre=idx_gre(:);
    else
        idx_gre=idx(sorted_ste.idx_shot(2,:));
    end
    ste_info=addvar2struct(ste_info,...
                           'idx_gre','idx','nc',...
                           'nv_cyc','nvf','nv',...
                           'resx_ste','resy_ste','resz_ste',...
                           'mot_cluster');
    if field_true(par,'cluster_b0')
        ste_info=addvar2struct(ste_info,...
                               'c1','c_b0_cent',...
                               'dist_max','dist_tot');
    end
    if pcaflag
        ste_info=addvar2struct(ste_info,...
                                   'scores',...
                                   'coeffs');
    end
end
