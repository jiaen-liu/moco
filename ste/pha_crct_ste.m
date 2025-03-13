function d=pha_crct_ste(d,pe,nstd,frac)
    pe=struct2double(pe);
    si = size(d);
    nr = si(1);
    nshot = si(end);
    necho = si(2);
    nch = si(end-1);
    nkyz = pe.hl(1);
    nblpl=total(all(reshape(pe.k,[2*necho,nkyz/necho])==0,1));
    nshot_blpl_cyc = nkyz/necho/nblpl;
    npair = floor(necho/2);
    ref_pha = zeros(nr, npair);
    for ishot = 1:nshot
        ishot_cyc = mod(ishot-1, nshot_blpl_cyc)+1;
        if ishot_cyc==1
            % look for the blipless data in this shot
            icyc = floor((ishot-1)/nshot_blpl_cyc)+1;
            % always the first shot is blipless!!!
            dblpls = d(:,:,:,(icyc-1)*nshot_blpl_cyc+1);
            % calculate the odd-even phase difference
            ref_pha=dpOddEven(dblpls,2,frac);
        end
        d(:,1:2:end,:,ishot)=d(:,1:2:end,:,ishot)./exp(1i*ref_pha/2);
        d(:,2:2:end,:,ishot)=d(:,2:2:end,:,ishot).*exp(1i*ref_pha/2);
    end
end
