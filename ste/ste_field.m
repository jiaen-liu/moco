function [ste_info,db0_fit]=ste_field(ste_info,db0,para,sorted_ste,par,mask)
% cluster the motion parmeters into a few groups
% within each, there is limited motion and the field
% changes can be approximated up to 1st order changes.
% create markers for motion corrupted full aquisitions
    fcrptd=ste_info.fcrptd;
    nv_cyc=ste_info.nv_cyc;
    % mot_par_f=ste_info.mot_par_f;
    nvf=ste_info.nvf;
    [nx,ny,nz,~,~]=size(db0);
    nc=ste_info.nc;
    ord=par.ord;
    gb0=par.gb0;
    
    if isfield(par,'pcaflag')
        pcaflag = par.pcaflag; 
        npca = par.npca; %To be provided, maybe add automatic choice of number of PCA Modes Later
    else
        pcaflag = 0;
    end
    %
    % save(rp('db0_comb.mat'),'db0_comb');
    % prepare the data struction for gre reconstruction
    % total number of shots in the gre data
    s2=para.sense_rate_s;
    s1=para.sense_rate_p;
    % nshot=para.n_partitions/s2*para.n_interleaves;
    nshot=sorted_ste.nshot;
    idx_shot_cenpf=sorted_ste.idx_shot_cenpf;
    nv=numel(idx_shot_cenpf);
    idx_gre=ste_info.idx_gre;
    cshot=ste_info.cshot;
    % resolution of ste
    resx_ste=ste_info.resx_ste;
    resy_ste=ste_info.resy_ste;
    resz_ste=ste_info.resz_ste;
    x_ste=([1:nx]-(1+nx)/2)*resx_ste;
    y_ste=([1:ny]-(1+ny)/2)*resy_ste;
    z_ste=([1:nz]-(1+nz)/2)*resz_ste;
    % image coordinate
    [xx_ste,yy_ste,zz_ste]=ndgrid(x_ste,y_ste,z_ste);
    xx_ste=xx_ste(mask);
    yy_ste=yy_ste(mask);
    zz_ste=zz_ste(mask);
    if ~field_true(par,'cluster_b0')
        % fit spherical harmonics to the field changes for each
        % accelerated image
        c_db0=sphere_harm_model_3d(db0(:,:,:,1:nv),x_ste,y_ste,z_ste,ord,mask);
        %
        if gb0
            A1=gen_spher_harm_poly(xx_ste,yy_ste,zz_ste,1);
        else
            A1=gen_spher_harm_poly(xx_ste,yy_ste,zz_ste,0);
        end
        A=gen_spher_harm_poly(xx_ste,yy_ste,zz_ste,ord);

        ste_info=setfield(ste_info,'c_db0',c_db0(:,ste_info.idx));
        % interplate the coefficients to each shot
        if ~sorted_ste.combined
            c_db0_gre=interp1(idx_shot_cenpf,c_db0.',[1:nshot],'pchip').';
        else
            c_db0_gre=c_db0(:,sorted_ste.idx_shot(2,:));
        end
        md=cell(nc,1);
    end

    if field_true(par,'cluster_b0')
% $$$         if isfield(par,'regress_nav') && par.regress_nav>0
% $$$             tmp=ste_info.c1;
% $$$             tmp=tmp.';
% $$$             nnav=size(tmp,1);
% $$$             if isfield(para,'loop_order') && ...
% $$$                     strcmp(para.loop_order,'zy_order') && ...
% $$$                     para.isgre
% $$$                 cyc_regress=para.n_partitions/...
% $$$                     para.sense_rate_s;
% $$$             else
% $$$                 cyc_regress=para.n_interleaves;
% $$$             end
% $$$             for i=1:size(tmp,2)
% $$$                 tmp(:,i)=regress_harm(tmp(:,i),...
% $$$                                           nnav*(sorted_ste.nshot_v+1),...
% $$$                                           cyc_regress,...
% $$$                                           par.regress_nav,...
% $$$                                           sorted_ste.nshot_v+1);
% $$$             end
% $$$             tmp=tmp.';
% $$$             ste_info.c1=tmp;
% $$$         end
        if isfield(par,'rect_sliding_window') && par.rect_sliding_window>0
            tmp=ste_info.c1;
            [n1,n2]=size(tmp);
            tmp2=tmp(:,1:floor(n2/nv_cyc)*nv_cyc);
            tmp2=reshape(tmp2,[n1,nv_cyc,floor(n2/nv_cyc)]);
            for i=1:n1
                tmp2(i,:,:)=ste_rect_sliding_window(reshape(tmp2(i,:,:),[nv_cyc,floor(n2/nv_cyc)]),...
                                                    par.rect_sliding_window);
            end
            tmp2=reshape(tmp2,[n1,nv_cyc*floor(n2/nv_cyc)]);
            tmp(:,1:nv_cyc*floor(n2/nv_cyc))=tmp2;
            if size(tmp,2)>nv_cyc*floor(n2/nv_cyc)
                tmp(:,nv_cyc*floor(n2/nv_cyc)+1:end)=repmat(tmp2(:,end),...
                                                            [1,size(tmp,2)-nv_cyc*floor(n2/nv_cyc)]);
            end
            c1=tmp;
            c1=c1(:,sorted_ste.idx_shot(2,:));
        else
            c1=ste_info.c1(:,sorted_ste.idx_shot(2,:));
        end
    end
    %Interpolate PCA coefficients and modes
    if pcaflag
       coeffs = ste_info.coeffs;
       scores = ste_info.scores;
       [nxs,nys,nzs,~,nte] = size(scores);
       nscores = nxs*nys*nzs;
       %nps = ny*nz;
       coefsinterp = zeros(nshot,npca,nte);
       scorefilt = reshape(scores,[nscores,npca,nte]);
       for echo=1:nte
           for mode = 1:npca
               %conjugate of coefficient is taken to account for the fact
               %b0=score*coef' which includes conjugating coef. but the
               %reconstructor will.* the fourier transform by coefsinterp thus
               %losing the conjugation.
               coefsinterp(:,mode,echo) = interp1(idx_shot_cenpf,conj(coeffs(:,mode,echo)),[1:nshot],'pchip').';
               %Interpolation of scores is likely also necessary. DOUBLE CHECK
           end
       end
       ste_info.coefsinterp = coefsinterp;
       ste_info.scorefilt = scorefilt;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:nc
        idx_gre_cl=find(idx_gre==i);
        if field_true(par,'cluster_b0')
            md{i}.db0=c1(1,idx_gre_cl);
            md{i}.gb0=c1(2:4,idx_gre_cl)*1e3;
            md{i}.c_db0=ste_info.c_b0_cent(:,i);
            md{i}.b0=[];
            if ~gb0
                md{i}.gb0(:)=0;
            end
        else
            % for each cluster, calculate the field difference to the
            % cluster reference up to the first order
            if gb0
                c_dif=zeros(4,numel(idx_gre_cl));
            else
                c_dif=zeros(1,numel(idx_gre_cl));
            end
            for j=1:numel(idx_gre_cl)
                c_dif(:,j)=A1\(A*(c_db0_gre(:,idx_gre_cl(j))-...
                                  c_db0_gre(:,cshot(i))));
            end
            if gb0
                % fit the difference of B0 as a gradient field within
                % each motion group
                tmp=c_dif;
                c_dif(1,:)=c_dif(1,:)*(0.5*(1/pi)^0.5);
                c_dif(2,:)=c_dif(4,:)*(3/4/pi)^0.5*1e3;
                c_dif(3:4,:)=tmp(2:3,:)*(3/4/pi)^0.5*1e3;
                md{i}.gb0=c_dif(2:4,:);
            else
                c_dif=c_dif*(0.5*(1/pi)^0.5);
                md{i}.gb0=zeros(3,numel(idx_gre_cl));
            end
            md{i}.db0=c_dif(1,:);
            md{i}.c_db0=c_db0_gre(:,cshot(i));
            md{i}.b0=db0(:,:,:,sorted_ste.idx_shot(2,cshot(i)));
        end
    end
    for i=1:nc
        tmp=ste_info.md{i};
        gb0_tmp=md{i}.gb0;
        db0_tmp=md{i}.db0;
        c_tmp=md{i}.c_db0;
        b0_tmp=md{i}.b0;
        tmp=setfield(tmp,'gb0',gb0_tmp);
        tmp=setfield(tmp,'db0',db0_tmp);
        tmp=setfield(tmp,'c_db0',c_tmp);
        tmp=setfield(tmp,'b0',b0_tmp);
        ste_info.md{i}=tmp;
    end
end
