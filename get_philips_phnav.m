function [data,dp1d,dp_fit,df_fit,df0,c,para]=get_philips_phnav(mid)
    frac=0.3;
    ord=2;
    [knav,label]=read_raw_philips(mid,'type',7);
    para=extract_para_philips(mid);
    nch=para.n_channels;
    nnav=para.nav1d_enable;
    nx=size(knav,1);
    nshot=numel(knav)/nx/nch/nnav;
    ishot_ref=floor(nshot/2);
    knav=reshape(knav,[nx,nch,nnav,nshot]);
    % fov correction
    lin_pha_ro=([0:nx-1].'-floor(nx/2))*2*pi*para.m_shift(1)/para.fovr/2;
    knav=knav.*exp(-1i*lin_pha_ro);
    data=fftmr(knav,-1,1);
    data=data(idx_truncate(nx,nx/2),:,:,:);
    nx=nx/2;
    dp0=angle(squeeze(sum(sum(data.*conj(data(:,:,:,ishot_ref)),2),1)));
    dp0=dp0(:);
    dp0uw=unwrap(dp0);
    dp0=dp0uw-dp0uw(ishot_ref)+dp0(ishot_ref);
    dp1d=angle(squeeze(sum(data.*conj(data(:,:,:,ishot_ref)),2)));
    dp1d=reshape(dp1d,[nx,nnav,nshot]);
    dp1duw=unwrap(dp1d,[],3);
    dp1duw=dp1duw-dp1duw(:,:,ishot_ref)+dp1d(:,:,ishot_ref);
    rms=sum(abs(data(:,:,end,ishot_ref)).^2,2).^0.5;
    mask=mask1d(rms,frac);
        idx_mask=find(mask);
    span_mask=idx_mask(end)-idx_mask(1);
    rel=0.5;
    while true
        maskDel=mask1d(rms,frac,rel);
        idx_mask_del=find(maskDel);
        span_mask_del=idx_mask_del(end)-idx_mask_del(1);
        if span_mask_del/span_mask>0.4
            break;
        else
            rel=rel-0.05;
        end
    end
    dp_fit=zeros(nx,nnav,nshot);
    c=zeros(ord+1,nnav,nshot);
    for i=1:nshot
        for j=1:nnav
            [dp_fit(:,j,i),c(:,j,i)]=polypha1d(dp1d(:,j,i),mask,ord,maskDel);
        end
    end
    df_fit=sum(dp_fit./2/pi.*reshape(para.te_nav1d*1e-3,[1,nnav]),2)/...
           sum((para.te_nav1d(:)*1e-3).^2,1);
    df_fit=reshape(df_fit,[nx,nshot]);
    df0=dp0/2/pi/para.te_nav1d/1e-3;
end

% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ clear;
% $$$ % $$$ cd ~/data_common/20240209_1;
% $$$ % $$$ mid=12;
% $$$ cd ~/data_common/20240227_1;
% $$$ mid=3;
% $$$ para=extract_para_philips(mid);
% $$$ % $$$ r=MRecon(mid2filename_philips(mid));
% $$$ % $$$ r.ReadData(7);
% $$$ % $$$ r.DcOffsetCorrection;
% $$$ % $$$ r.PDACorrection;
% $$$ % $$$ r.RandomPhaseCorrection;
% $$$ % $$$ r.MeasPhaseCorrection;
% $$$ [knav,label]=read_raw_philips(mid,'type',7);
% $$$ 
% $$$ nch=para.n_channels;
% $$$ nshot=para.n_main_tr;
% $$$ nnav=para.nav1d_enable;
% $$$ nte=length(para.te_contr);
% $$$ 
% $$$ [kd,para,l,cov_mat]=recon_bmir_epi(mid,'k_return',1);
% $$$ 
% $$$ kc=kd(:,:,:,:,find(l.ky(1:nte:end)==0&l.kz(1:nte:end)==0));
% $$$ prof=fftmr(kc,-1,1);
% $$$ 
% $$$ 
% $$$ knav=reshape(knav,[nx,nch,nnav,nshot]);
% $$$ 
% $$$ nx=size(knav,1);
% $$$ lin_pha_ro=([0:nx-1].'-floor(nx/2))*2*pi*para.m_shift(1)/para.fovr/2;
% $$$ knav=knav.*exp(-1i*lin_pha_ro);
% $$$ data=fftmr(knav,-1,1);
% $$$ data=data(idx_truncate(nx,nx/2),:,:,:);
% $$$ nx=nx/2;
% $$$ % correct fov
% $$$ 
% $$$ 
% $$$ dp=angle(sum(squeeze(sum(data.*conj(data(:,:,:,floor(nshot/2))),2)),1));
% $$$ dp1d=angle(squeeze(sum(data.*conj(data(:,:,:,floor(nshot/2))),2)));
% $$$ dp=squeeze(dp).';
% $$$ dp_ori=dp;
% $$$ for i=1:nnav
% $$$     dp(:,i)=regress_harm(dp_ori(:,i),...
% $$$                          nshot,para.n_interleaves,...
% $$$                          20,1);
% $$$ end
% $$$ 
% $$$ df=dp./2/pi./para.te_nav1d(:).'/1e-3;
% $$$ 
% $$$ data_sum=reshape((sum(abs(data).^2,2).^0.5),[nx,nnav,nshot]);
% $$$ 
% $$$ for i=1:nnav
% $$$     figure;plot(sz(data_sum(:,i,:)));setLineColor;
% $$$ end
% $$$ figure;plot(df)
% $$$ 
% $$$ % $$$ nt=size(dp,1);
% $$$ % $$$ A=[eye(nt,nt)*para.te_nav1d(1)*2*pi*1e-3;...
% $$$ % $$$    eye(nt,nt)*para.te_nav1d(2)*2*pi*1e-3];
% $$$ % $$$ df_fit=A\[dp(:,1);dp(:,2)];
% $$$ % $$$ figure;plot(df_fit);ylim([-0.5,4.5])
% $$$ % $$$ dp_clean=dp-df_fit*2*pi.*para.te_nav1d(:).'*1e-3;