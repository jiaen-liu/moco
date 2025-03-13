function mkps=sense_mask_k(np,ns,sp,ss,ds,dy,dz)
    mkps=false(np,ns);
    %% shift kz
    mks=false(ns,1);
    mks(1:ss:ns)=true;
    for ip=1:sp:np
        mkps(ip,:)=mks;
        mks=circshift(mks,ds);
    end
    mkps=circshift(mkps,dy,1);
    mkps=circshift(mkps,dz,2);
end