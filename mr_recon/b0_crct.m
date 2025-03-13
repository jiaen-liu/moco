function d=b0_crct(d,nav,ste,te,para)
    use_ste=0;
    if isempty(nav) & ~isempty(ste)
        use_ste=1;
    elseif isempty(nav) & isempty(ste)
        warning('*** Need navigator or STE data for B0 correction! ***');
        return;
    end
    ntr=para.n_main_tr;
    nintl=para.n_interleaves;
    if use_ste
        freq=ste.df_shot;
        freq=reshape(freq,[1,1,1,1,nintl,ntr/nintl]);
    else
        dp=angle(sum(sum(nav.d.*...
                         conj(nav.d(:,1,:,:,floor(ntr/2))),...
                         4),...
                     1));
        freq=dp/2/pi/nav.te/1e-3;
        freq=reshape(freq,[1,1,para.n_slices,...
                           1,nintl,ntr/nintl]);
    end
    si=size(d);
    d=reshape(d,[si(1:end-1),nintl,ntr/nintl]);
    te=reshape(te,[1,length(te)]);
    if para.int_te_shift
        dte=para.echo_spacing*1e-3/nintl;
        te=te+reshape([0:nintl-1]*dte,[1,1,1,1,nintl]);
    end
    d=d./exp(1i*2*pi*1e-3*(freq.*te));
end