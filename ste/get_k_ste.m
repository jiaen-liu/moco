function dk=get_k_ste(sorted_ste,iv,nshot_v)

    [nr,necho,nch,nshot]=size(sorted_ste.d_nblpl);
    np=sorted_ste.dims_k(2);
    ns=sorted_ste.dims_k(3);
    ncontr=sorted_ste.ncontr;
    necho_contr=necho/ncontr;
    nv=length(iv);
    dk=zeros(nr,np,ns,nch,ncontr,nv);
    for ivv=1:nv
        ishot=(iv(ivv)-1)*nshot_v+[1:nshot_v];
        ishot_k=mod(ishot-1,size(sorted_ste.k_nblpl,3))+1;
        % ishot_k=(ivp-1)*nshot_v+[1:nshot_v];
        k=sorted_ste.k_nblpl(:,:,ishot_k);
        dk_raw=sorted_ste.d_nblpl(:,:,:,ishot);
        for icontr=1:ncontr
            for ie=1:necho_contr
                for is=1:nshot_v
                    dk(:,floor(np/2)+1+k(1,ie,is),...
                       mod(floor(ns/2)+k(2,ie,is),ns)+1,:,icontr,ivv)=...
                        dk_raw(:,ie+(icontr-1)*necho_contr,:,is);
                end
            end
        end
    end
end