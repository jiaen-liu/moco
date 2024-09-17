function y=fov_crct(data,k,para)
    if nargin<4
        is=1;
    end
    sign_ph = -1.0;
    sign_sl = -1.0;
    [nx,n_echo_shot,ns,nch,ntr]=size(data);
    y=zeros(size(data));
    % [npe,n_echo_k,ntrk]=size(k);
    npe=size(k,1);
    k=reshape(k,[npe,n_echo_shot,ns,1,ntr]);
% $$$     if n_echo_k~=n_echo_shot || ntr~=ntrk
% $$$         error('*** Array dimensions do not match! ***');
% $$$     end
    for is=1:ns
        coor_new = gen_coordinate(para.snormal(:, is),...
                                  para.prot);
        % correct phase encoding direction
        p_shift = [para.x_shift(is),...
                   para.y_shift(is),...
                   para.z_shift(is)]*...
                  coor_new(:, 2);
        delp = p_shift/para.fovp;
        % correct slice direction
        s_shift = [para.x_shift(is),...
                   para.y_shift(is),...
                   para.z_shift(is)]*...
                  coor_new(:, 3);
        dels = s_shift/(para.sthickness*...
                        (1+para.sl_oversamp));
        y(:,:,is,:,:)=data(:,:,is,:,:).*...
            exp(sign_ph*2*pi*1i*delp*k(1,:,is,:,:)).*...
            exp(sign_sl*2*pi*1i*dels*k(2,:,is,:,:));
    end
end