function b1info=gen_b1_info(mid,cov_mat)
    para=extract_para(mid);
    coord=get_coordinate(para,1,0);
    b1data=recon_amri_epi(mid,'k_return',1);
    b1data=squeeze(b1data(:,:,:,:,:,1));
    [nx,ny,nz,nch]=size(b1data);
    if para.dimen==3
        b1data=fftmr(b1data,-1,[1,2,3]);
    else
        b1data=fftmr(b1data,-1,[1,2]);
    end
    mask=auto_im_mask(sum(abs(b1data),4));
    [sensit,ref]=sense_m(b1data,mask,cov_mat);
    b1info.mask=mask;
    b1info.sensit=sensit;
    b1info.ref=ref;
    b1info.coord=coord;
    fn_b1=fullfile(pwd,['mid',num2str(mid),'.b1info.mat']);
    save(fn_b1,'b1info');
end