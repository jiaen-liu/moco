function db0_comb=ste_comb_b0(db0f,db0,ste_info,par,sorted_ste)
    nv_cyc=sorted_ste.nparv_cyc;
    si=size(db0);
    nx=si(1);
    ny=si(2);
    nz=si(3);
    nv=si(end);
    if ~isfield(par,'mot_cmb') || par.mot_cmb==0
        % use fast navigator
        db0_comb=db0;
    elseif par.mot_cmb==2
        % use full navigator
        db0_comb=repmat(reshape(db0f,[nx,ny,nz,1,nv]),[1,1,1,nv_cyc,1]);
    elseif par.mot_cmb==1
        % combine
        db0_comb=db0;
        db0_comb(:,:,:,:,~ste_info.fcrptd)=repmat(...
            reshape(db0f(:,:,:,~ste_info.fcrptd),[nx,ny,nz,1,total(~ste_info.fcrptd)]),...
            [1,1,1,nv_cyc,1]);
    end
    % the first volume is not stable
    db0_comb(:,:,:,:,1)=db0_comb(:,:,:,:,2);
end
