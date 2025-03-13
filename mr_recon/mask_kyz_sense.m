function mkps=mask_kyz_sense(ny,nz,s1,s2,caipi,dy,dz)
mkps=false(ny,nz);
%% shift kz
mks=false(nz,1);
mks(1:s2:nz)=true;
for iy=1:s1:ny
    mkps(iy,:)=circshift(mks,mod((iy-1)/s1*caipi,s2));
end
mkps=circshift(mkps,dy,1);
mkps=circshift(mkps,dz,2);
end