function dste=detrend_ste(dste, ishot_ref, par_regress)
% first order correction of frequency drift
% maybe zero order is better when motion happens
    if nargin<3
        par_regress.enable=0;
    end
    d = fftmr(dste.d,-1,1);
    si = size(d);
    nr = si(1);
    x=[0:nr-1].';
    nshot = si(end);
    nch = si(end-1);
    necho = si(2);
    nkyz = dste.pe.hl(1);
    nshot_kyz = double(nkyz)/necho;
    kyz = reshape(dste.pe.k, [2, necho, nshot_kyz]);
    blplsm = col(sum(sum(abs(kyz), 1), 2)) == 0;
    nshot_cyc = round(nshot_kyz/total(blplsm));
    % get the blipless data of the ref shot
    if nargin<2 
        ishot_ref = floor(nshot/2)+1;
    end
    icyc_ref = floor((ishot_ref-1)/nshot_cyc)+1;
    % always the first shot is blipless!!!
    dblpls_ref = d(:,:,:,(icyc_ref-1)*nshot_cyc+1);
    mag = sum(abs(dblpls_ref(:,1,:)).^2, 3).^0.5;
    
    % generate a mask based on the magnitude
    mask = mask1d(mag,0.3);
    xm = x(find(mask));
    nm = numel(xm);
    te = dste.te(:)*1e-3;
    te_odd = te(1:2:end);
    te_even = te(2:2:end);
    dte_odd = te_odd(2)-te_odd(1);
    nblpl = ceil(nshot/nshot_cyc);
    df_blpl = zeros(nblpl, 2);
    iblpl = 1;
    pi2 = 2*pi;
    % calculate frequency of reference
    nodd = ceil(necho/2.0);
    nhecho = nodd;
    dblpls_ref_odd = dblpls_ref(:,1:2:end,:);
    dblpls_ref_even = dblpls_ref(:,2:2:end,:);
    dp_ref_odd=angle(sum(sum(dblpls_ref_odd.*conj(dblpls_ref_odd(:,1,:)),1),3));
    dp_ref_odd=dp_ref_odd(:);
    fref=fit_freq_linear(dp_ref_odd,te_odd);
    
    if par_regress.enable
        mask_te=true(length(te_odd),1);
    else
        mask_te=te_odd>5e-3&te_odd<10e-3;
    end
    
    if all(~mask_te)
        mask_te(end) = true;
    end
    
    dblplsm_all=false(nshot,1);
    dblplsm_all(1:nshot_cyc:end)=true;
    dblpls_all=d(:,:,:,find(dblplsm_all));
    dp=squeeze(angle(sum(sum(dblpls_all.*conj(dblpls_ref),1),3)));
    if par_regress.enable
        for i=1:necho
            if isfield(dste.para,'loop_order') && ...
                    strcmp(dste.para.loop_order,'zy_order') && ...
                    dste.para.isgre
                dp(i,:)=regress_harm(col(dp(i,:)),...
                                         nshot,dste.para.n_partitions/dste.para.sense_rate_s,...
                                         par_regress.order,nshot_cyc);
            else
                dp(i,:)=regress_harm(col(dp(i,:)),...
                                         nshot,dste.para.n_interleaves,...
                                         par_regress.order,nshot_cyc);
            end

        end
    end
    dp_odd=dp(1:2:end,:);
    df_odd=col(te_odd(find(mask_te))*pi2)\...
                   dp_odd(find(mask_te),:);
    dp_even=dp(2:2:end,:);
    df_even=col(te_even(find(mask_te))*pi2)\...
                   dp_even(find(mask_te),:);
    df_blpl(:,2)=0.5*(df_odd(:)+df_even(:));
    df_blpl(:,1)=find(dblplsm_all);
% $$$     for ishot = 1:nshot
% $$$         ishot_cyc = mod(ishot-1, nshot_cyc)+1;
% $$$         if ishot_cyc==1
% $$$             % look for the blipless data in this shot     
% $$$             icyc = floor((ishot-1)/nshot_cyc)+1;
% $$$             % always the first shot is blipless!!!
% $$$             dblpls = d(:,:,:,(icyc-1)*nshot_cyc+1);
% $$$             % zero order
% $$$             dblpls_odd = dblpls(:, 1:2:end, :);
% $$$             dblpls_even = dblpls(:, 2:2:end, :);
% $$$             
% $$$             dp_odd = angle(sum(sum(dblpls_odd.*conj(dblpls_ref_odd), 3),1));
% $$$             dp_odd=dp_odd(:);
% $$$             dp_even = angle(sum(sum(dblpls_even.*conj(dblpls_ref_even), 3),1));
% $$$             dp_even=dp_even(:);
% $$$ 
% $$$             df_odd=col(te_odd(find(mask_te))*pi2)\...
% $$$                    col(dp_odd(find(mask_te)));
% $$$             df_even=col(te_even(find(mask_te))*pi2)\...
% $$$                    col(dp_even(find(mask_te)));
% $$$             df_blpl(iblpl,2)=0.5*(df_odd+df_even);
% $$$ % $$$             
% $$$ % $$$             df_odd=dp_odd./te_odd/pi2;
% $$$ % $$$             df_even=dp_even./te_even/pi2;
% $$$ % $$$             df_blpl(iblpl, 2)=mean(0.5*(df_odd(find(mask_te))+...
% $$$ % $$$                                         df_even(find(mask_te))));
% $$$             df_blpl(iblpl, 1)=ishot-1;
% $$$             iblpl=iblpl+1;
% $$$         end
% $$$     end
% $$$     if par_regress.enable
% $$$         df_blpl(:,2)=regress_harm(df_blpl(:,2),...
% $$$                                       nshot,dste.para.n_interleaves,...
% $$$                                       par_regress.order,nshot_cyc);
% $$$     end
    df_blpl(:, 2)=df_blpl(:, 2)+fref;
    % remove field fluctuation from each shot
    df_shot = interp1(df_blpl(:,1), df_blpl(:,2), ...
                      [0:nshot-1].','spline');
    d=d./exp(1i*pi2*(reshape(te,[1,necho]).*reshape(df_shot,[1,1,1,nshot])));
    dste.d = fftmr(d,1,1);
    dste.df_shot=df_shot;
    dste.df_blpl=df_blpl;
end

