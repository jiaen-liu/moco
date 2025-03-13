function [d,te]=prep_nav(mid,navtype,para,no_regress)
    if nargin<4
        no_regress=0;
    end
    d=readSortSiem(mid,'nav',1,'main',1);
    [dmdh,mdhTempl]=readmdh(mid,'noise');
    if navtype==0
        return;
    end
    i_nav = conditional(navtype==0,1,navtype);
    navmdh=cast2struct(dmdh(:,1,i_nav),mdhTempl);
    % determine the segment index of the navigator data
    navmdh_cur=cast2struct(dmdh(:,1,1),mdhTempl);
    idx1=1;
    idx2=navmdh_cur.samples;
    for i=1:i_nav-1
        navmdh_cur=cast2struct(dmdh(:,1,i+1),mdhTempl);
        idx1=idx2+2;
        idx2=idx2+navmdh_cur.samples;
    end
    d=d(idx1:idx2,1,:,:,:);
    % correct readout polarity
    if ~evalmaskbit(navmdh,25)
        d=flipdim(d,1);
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
end