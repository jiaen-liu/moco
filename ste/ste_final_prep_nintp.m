% 20230621: Jiaen Liu, normalize b1 data by the covariance matrix to improve sensitivity estimation, especially useful for 10.5T
function [para,ste_info,b0_fit,mot_par_test,b1]=ste_final_prep_nintp(sorted_ste,ste_info,mask,mask_brain,db0,mot_f,...
                                                      imf_uncomb,par)
% prepare the parameters, sense and b0 maps to be used in reconstructing gre images
    para=sorted_ste.para;
    para.idx_shot_cenpf=sorted_ste.idx_shot_cenpf;
    % combine the field changes from full-fov and
    % accelerated images
    % Group motion parameters to clusters based on k-means
    ste_info=ste_field(ste_info,db0,...
                       para,sorted_ste,par,mask_brain);
    if par.discard
        ste_info=ste_discard(ste_info,sorted_ste);
    end
    % coordinate of the ste image
    nx_ste=para.steref_dim_r;
    ny_ste=para.steref_dim_p;
    nz_ste=para.steref_dim_s;
    res_ste=[para.steref_res_r,...
             para.steref_res_p,...
             para.steref_res_s];
    x_ste=([1:nx_ste]-(1+nx_ste)/2)*res_ste(1);
    y_ste=([1:ny_ste]-(1+ny_ste)/2)*res_ste(2);
    z_ste=([1:nz_ste]-(1+nz_ste)/2)*res_ste(3);
    [x_ste,y_ste,z_ste]=ndgrid(x_ste,y_ste,z_ste);
    coord_ste=cat(4,x_ste,y_ste,z_ste);
    % coordiante of the main acquision
    nx=para.nr;
    ny=para.np;
    nz=para.n_slices*para.n_partitions;
    x=([1:nx]-(nx+1)/2)*para.resr;
    y=([1:ny]-(ny+1)/2)*para.resp;
    z=([1:nz]-(nz+1)/2)*para.ress;
    [x,y,z]=ndgrid(x,y,z);
    coord_main=cat(4,x,y,z);
    % calculate interpolation kernel for B0
    % which is always based on ste
    intp_ker_ste=mkl_interp_kernel(coord_ste,coord_main);
    % interpolation kernel for B1 depends on the use of
    % external reference
    intp_ker_ext=[];
    % ----------------------------------------------- %
    % Obtain sensitivity maps for each cluster
    % ----------------------------------------------- %
    isen=zeros(ste_info.nc,1);
    for i=1:ste_info.nc
        isen(i)=ste_info.md{i}.isen;
    end
    if field_true(par,'use_pri_b1')
        b1=interp3_nmat(coord_main,coord_ste,read_data(rp('sensit_pri.svd')));
    else
        if isfield(par,'use_ext_ref') && ...
                ~isempty(par.use_ext_ref) && ...
                par.use_ext_ref(1)>0
            % use sense reference from a seperate scan
            % par.use_ext_ref(1) is a mid number
            % par.use_ext_ref(2) is mask_threshold
            b1=repmat(ste_sense_ext(par,sorted_ste),[1,1,1,1,ste_info.nc]);
            % calculate interpolation kernel for B1 
            % based on the external reference
            mid_pimg=par.use_ext_ref(1);
            para_pimg=getfield(sorted_ste.para_pimg,['mid' num2str(mid_pimg)]);
            coord_ext_magnet=get_coordinate(para_pimg,1,0,0);
            coord_main_magnet=get_coordinate(para,1,0,0);
            intp_ker_ext=mkl_interp_kernel(coord_ext_magnet,coord_main_magnet);
        else
            % use internal reference scan
            if isempty(imf_uncomb) || ...
                    isempty(mot_f)
                error('*** No image data is provided for internal reference of B1! ***');
            end
            b1=ste_sense(imf_uncomb,mask,...
                         [sorted_ste.te(1+floor(end/4)),...
                          sorted_ste.te(1+end/2+floor(end/4))]*1e-3,...
                         mot_f(:,isen),[],[],par,res_ste/par.ste_intp_res);
        end
    end
    mask_b1=squeeze(sum(abs(b1),4))>0;
    % need to interplate mask_b1 to ste
    if isfield(par,'use_ext_ref') && ...
                ~isempty(par.use_ext_ref) && ...
                par.use_ext_ref(1)>0
        coord_ste_magnet=get_coordinate(para,1,1,0);
        mask_b1=interp3_nmat(coord_ext_magnet,coord_ste_magnet,single(mask_b1))>0.5;
    end
    nch=para.n_channels;
    [nxb1,nyb1,nzb1,~,~]=size(b1);
    nb1=numel(b1)/nch/ste_info.nc;
    b1=reshape(b1,[nb1,nch,ste_info.nc]);
    b1=permute(b1,[1,3,2]);
    b1=reshape(b1,[nb1*ste_info.nc,nch]);
    
    if ~(isfield(par,'use_ext_ref') && ...
         ~isempty(par.use_ext_ref) && ...
         par.use_ext_ref(1)>0)
        % for external sense reference scan, this was already done and not needed here.
        % see ste_sense_ext.m for details
        b1=b1*conj(chol(sorted_ste.inv_cov,'lower'));
    end
    b1=reshape(b1,[nb1,ste_info.nc,nch]);
    b1=permute(b1,[1,3,2]);
    b1=reshape(b1,[nxb1,nyb1,nzb1,nch,ste_info.nc]);
    b1=b1/prctile(abs(b1(:)),95);
    para.b1n=b1;
    % ----------------------------------------------- %
    % B0 maps for each cluster
    % ----------------------------------------------- %
    b0=zeros(nx_ste,ny_ste,nz_ste,ste_info.nc);
    mask_ext=zeros(nx_ste,ny_ste,nz_ste);
    mask_brain_ext=zeros(nx_ste,ny_ste,nz_ste);
    % shrink mask_brain
    mask_brain=volerode(mask_brain,2);
    kernel = fspecial('gaussian',7,1);
    for i = 1:nz_ste
        mask_ext(:, :, i) = imfilter(double(mask(:, :, i)), kernel);
        mask_brain_ext(:,:,i)=imfilter(double(mask_brain(:, :, i)), kernel);
    end
% $$$     mask_b1_ext=zeros(nx_ste,ny_ste,nz_ste,ste_info.nc);
% $$$     for j=1:size(mask_b1,4)
% $$$         for i=1:nz_ste
% $$$             mask_b1
% $$$         end
% $$$     end
    mask_ext(find(mask)) = 1.0;
    Am=gen_spher_harm_poly(x_ste,y_ste,z_ste,par.ord);

    
    for i=1:ste_info.nc
        b0(:,:,:,i)=reshape(Am*ste_info.md{i}.c_db0,...
                            [nx_ste,ny_ste,nz_ste]);
    end
    % b0=b0.*mask_ext;
    b0=b0.*mask_b1;
    
% $$$     b0_meas=zeros(nx_ste,ny_ste,nz_ste,ste_info.nc);
% $$$     for i=1:ste_info.nc
% $$$         b0_meas(:,:,:,i)=ste_info.md{i}.b0;
% $$$     end
    % combine measured and fitted b0
    % b0=b0.*(1-mask_brain_ext)+b0_meas.*mask_brain_ext;
    % b0=b0.*mask_main_ext;

    b0_fit=zeros(size(db0));
    b0_fit=combine_dim(b0_fit,[4,5]);
    nshot=sorted_ste.nshot;
    idx_gre=ste_info.idx_gre;
    idx=ste_info.idx(1:length(sorted_ste.idx_shot_cenpf));
    db0tmp=zeros(nshot,1);
    gb0tmp=zeros(3,nshot);
    mtmp=zeros(3,2,nshot);
    for i=1:ste_info.nc
        gb0tmp(:,idx_gre==i)=ste_info.md{i}.gb0;
        db0tmp(idx_gre==i)=ste_info.md{i}.db0;
        mtmp(:,:,idx_gre==i)=ste_info.md{i}.m;
    end
    mot_par_test=zeros(6,length(sorted_ste.idx_shot_cenpf));
    A1=[x_ste(:),y_ste(:),z_ste(:)]*1e-3;
    for i=1:ste_info.nc
        sc=total(idx==i);
        b0_fit(:,:,:,idx==i)=...
            b0(:,:,:,i)+...
            reshape(db0tmp(sorted_ste.idx_shot_cenpf(idx==i)),...
                    [1,1,1,sc])+...
            reshape(A1*...
                    gb0tmp(:,sorted_ste.idx_shot_cenpf(idx==i)),...
                    [nx_ste,ny_ste,nz_ste,sc]);
        mot_par_test(:,idx==i)=...
            reshape(mtmp(:,:,sorted_ste.idx_shot_cenpf(idx==i)),...
                    [6,sc]);
    end
    if field_true(par,'use_pri_b0')
        para.b0=interp3_nmat(coord_main,coord_ste,read_data(rp('b0_pri.svd')));
    else
        para.b0=b0;
    end
    
    para.j=5;
    para.nm=ste_info.nc;
    para.n_iter_reset=par.n_iter_reset;
    para.n_iter=par.n_iter;
    for i=1:ste_info.nc
        md(i)=ste_info.md{i};
    end
    para=setfield(para,'md',md);
    para.kos=par.kos;
    para.en_parfor=par.en_parfor;
    para.intp_ker_ste=intp_ker_ste;
    para.intp_ker_ext=intp_ker_ext;
end