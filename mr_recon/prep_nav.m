function [d,te,fov,res]=prep_nav(mid,navtype,para,no_regress,whitening)
    if nargin<4
        no_regress=0;
    end
    if nargin<5
        whitening=0;
    end
    d=readSortSiem(mid,'nav',1,'main',1);
    [dmdh,mdhTempl]=readmdh(mid,'noise');
    if navtype==0
        return;
    end
    n_echo_ste=0;
    if para.b_ste_en
        n_echo_ste=para.n_echo_steref;
    end
    i_nav = conditional(navtype==0,1,navtype);
    navmdh=cast2struct(dmdh(:,1,i_nav+n_echo_ste),mdhTempl);
    % determine the segment index of the navigator data
    navmdh_cur=cast2struct(dmdh(:,1,1+n_echo_ste),mdhTempl);
    idx1=1;
    idx2=navmdh_cur.samples;
    for i=1:i_nav-1
        navmdh_cur=cast2struct(dmdh(:,1,i+1+n_echo_ste),mdhTempl);
        idx1=idx2+1;
        idx2=idx2+navmdh_cur.samples;
    end
    d=d(idx1:idx2,1,:,:,:);
    % correct readout polarity
    if ~evalmaskbit(navmdh,25) && navtype==1
        % the direction of nav is opposite to that of image in PE and Slice direction
        % so if navtype is 2 or 3, do not flip
        % only flip for readout direction, then it will be consistent with recon_amri_epi
        % and prep_ste
        % also check gen_coordinate "x (readout) gets inverted in Jiaen's recon convention"
        d=flipdim(d,1);
    end
    if whitening
        cov_mat=covSiem(mid);
        d=covNorm(d,cov_mat,4);
    end
    % regress out echo-shift eddy current
    if para.int_te_shift && ~no_regress
        d=steRegress(d,para.n_interleaves,...
                     floor(para.n_interleaves/2));
    end
    nr=para.nav_dim(navtype);
    % regridding
    % 'ramp_samp_frac',para.ramp_samp_frac,...
    if isfield(navmdh,'rsamptime') && navmdh.rsamptime>0
        para_rsamp=struct('nr',nr,...
                          't_dwell',1/para.nav_bw/2*1e9,...
                          'ramp_dur',double(navmdh.ramp),...
                          'ramp_samp_frac',double(navmdh.rsamptime)/double(navmdh.ramp),...
                          'nr_os',nr*2);
        d=regridding_arr(d,para_rsamp,0.25);
    else
        nr=para.nav_dim(navtype);
        d=fftmr(d,-1,1);
    end
    d=d(idx_truncate(nr*2,nr),:,:,:,:);
    te=navmdh.te;
    switch navtype
      case 1
        fov=para.fovr;
      case 2
        fov=para.fovp;
      case 3
        fov=para.sthickness*2;
    end
    res=fov/nr;
    % correct fov off center
    % TBD: the correct shift
    if navtype>1
        coor = gen_coordinate(para.snormal(:, 1), para.prot);
        shift = [para.x_shift(1), para.y_shift(1), para.z_shift(1)]*...
                coor(:, navtype);
        d=fftmr(fftmr(d,1,1).*exp(-1i*([0:nr-1].'-floor(nr/2))*2*pi*shift/fov),-1,1);
    end
end
