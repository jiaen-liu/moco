function [kyz,mkps]=gen_pe_sense(ny,nz,s1,s2,ds,dy,dz)
mkps=false(ny,nz);
mkps(1:s1:end,1:s2:end)=1;
mkps=circshift(mkps,dy,1);
mkps=circshift(mkps,dz,2);

ky=[-floor(ny/2):floor((ny-1)/2)];
kz=[-floor(nz/2):floor((nz-1)/2)];
[ky,kz]=ndgrid(ky,kz);
ky=ky(mkps(:)).';
kz=kz(mkps(:)).';

kz=reshape(kz,[ny/s1,nz/s2]);
for iy=1:ny/s1
    kz(iy,:)=kz(iy,:)+mod((iy-1)*ds,s2);
end
kyz=[ky;kz(:).'];