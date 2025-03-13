function data=fov_crct_ste(data, k, para, is)
    if nargin<4
        is=1;
    end
    sign_ph = -1.0;
    sign_sl = -1.0;
    coor_new = gen_coordinate(para.snormal(:, is), para.prot);
    % correct phase encoding direction
    p_shift = [para.x_shift(is), para.y_shift(is), para.z_shift(is)]*...
              coor_new(:, 2);
    delp = p_shift/para.fovp;
    % correct slice direction
    s_shift = [para.x_shift(is), para.y_shift(is), para.z_shift(is)]*...
              coor_new(:, 3);
    dels = s_shift/(para.sthickness*(1+para.sl_oversamp));

    si = size(data);
    necho = si(2);
    nshot = si(end);
    nk = size(k,2);
    for i = 1:nshot 
        for j = 1:necho
            ik = mod((j-1)+(i-1)*necho,nk)+1;
            data(:, j, :, i) = data(:, j, :, i)*...
                (exp(sign_ph*1i*2*pi*k(1, ik)*delp(1))*...
                 exp(sign_sl*1i*2*pi*k(2, ik)*dels(1)));
        end
    end
end
