function [k_pimg,im_pimg,para_pimg]=ste_nav_ref(par,sorted_ste)
    if ~isfield(sorted_ste,'k_pimg')
        error(['*** Parallel imaging ', ...
               'reference was ',...
               'not provided! ***']);
    end
    if ~isfield(par,'ste_nav_ref')
        error(['*** MID of parallel ',...
               'imaging reference ',...
               'was not provided! ***']);
    end
    mid_pimg=par.ste_nav_ref(1);
    % get data and parameters
    if ~isfield(sorted_ste.k_pimg,['mid' num2str(mid_pimg)])
        error(['*** Parallel imaging data for MID '...
               num2str(mid_pimg) ...
               ' was not provided! ***']);
    end
    k_pimg=getfield(sorted_ste.k_pimg,['mid' num2str(mid_pimg)]);
    sik=size(k_pimg);
    k_pimg=k_pimg(:,:,:,:,1,1,1,1);
    im_pimg=fftmr(k_pimg,-1,[1,2,3])*prod(sik(1:3));
    % parameters
    para_pimg=getfield(sorted_ste.para_pimg,['mid' num2str(mid_pimg)]);
    pix_shift=[para_pimg.nr,...
               para_pimg.np,...
               para_pimg.n_partitions];
    pix_shift=-mod(1+pix_shift,2)/2;
    if field_true(par,'shift_ste')
        im_pimg=shift(im_pimg,pix_shift);
    end
    % interpolate external reference to ste
    nx_ste=sorted_ste.para.steref_dim_r;
    ny_ste=sorted_ste.para.steref_dim_p;
    nz_ste=sorted_ste.para.steref_dim_s;
    res_ste=[sorted_ste.para.steref_res_r,...
             sorted_ste.para.steref_res_p,...
             sorted_ste.para.steref_res_s];
    x_ste=([1:nx_ste]-(1+nx_ste)/2)*res_ste(1);
    y_ste=([1:ny_ste]-(1+ny_ste)/2)*res_ste(2);
    z_ste=([1:nz_ste]-(1+nz_ste)/2)*res_ste(3);
    [x_ste,y_ste,z_ste]=ndgrid(x_ste,y_ste,z_ste);
    coord_ste=cat(4,x_ste,y_ste,z_ste);
    
    nx_ref=para_pimg.nr;
    ny_ref=para_pimg.np;
    nz_ref=para_pimg.n_partitions;
    res_ref=[para_pimg.resr,...
             para_pimg.resp,...
             para_pimg.ress];
    x_ref=([1:nx_ref]-(1+nx_ref)/2)*res_ref(1);
    y_ref=([1:ny_ref]-(1+ny_ref)/2)*res_ref(2);
    z_ref=([1:nz_ref]-(1+nz_ref)/2)*res_ref(3);
    [x_ref,y_ref,z_ref]=ndgrid(x_ref,y_ref,z_ref);
    coord_ref=cat(4,x_ref,y_ref,z_ref);
    im_pimg=interp3_nmat(coord_ref,coord_ste,im_pimg);
end