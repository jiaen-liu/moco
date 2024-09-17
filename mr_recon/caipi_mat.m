function [m,I,J,V]=caipi_mat(nx,ny,nz,ry,rz,zshift,dy,dz)
% m is a sparse matrix
    n=nx*ny*nz;
    nr=ry*rz;
    I=reshape(repmat([1:n],[nr,1]),[nr*n,1]);
    J=zeros(nr,n);
    V=zeros(nr,n);
    % minimum ky or kz
    kymin=-floor(ny/2)/ny*2*pi;
    kzmin=-floor(nz/2)/nz*2*pi;
    % steps of ky or kz in full fov
    dky=2*pi/ny;
    dkz=2*pi/nz;
    nys=ny/ry;
    nzs=nz/rz;
    nxy=nx*ny;
    % ----------------------------------------------- %
    % aliasing as a function of z and zshift
    % dy(voxels)=-ny/(ry*rz)*(zj-zi)/(nz/rz)*zshift
    % ----------------------------------------------- %
    % ----------------------------------------------- %

    % J is the index of alised y and z positions
    % generate the first row of J, the rest is shifted
    % incrementally by 1
    for i=1:n
        x=mod(i-1,nx)+1;
        y=mod(floor((i-1)/nx),ny)+1;
        z=floor((i-1)/nxy)+1;
        for iz=1:rz
            for iy=1:ry
                ty=mod((iy-1)*nys-(iz-1)*zshift*ny/nr+y-1,ny);
                tz=mod((iz-1)*nzs+z-1,nz);
                J(iy+(iz-1)*ry,i)=x+ty*nx+tz*nxy;
            end
        end
    end
% $$$     for k=2:n
% $$$         J(:,k)=mod(J(:,k-1),n)+1;
% $$$     end
    % calculate V
    ky_start=kymin+dy*2*pi/ny;
    kz_start=kzmin+dz*2*pi/nz;
    for iz=1:rz
        for iy=1:ry
            V(iy+(iz-1)*ry,1)=exp(-1i*(ky_start/nr*ny*zshift*(iz-1)))...
                *exp(1i*dz*2*pi/nz*(iz-1)*nzs)...
                *exp(1i*dy*2*pi/ny*(iy-1)*nys)...
                *exp(-1i*2*pi*floor(ny/2)*(iy-1)/ry)...
                *exp(-1i*2*pi*floor(nz/2)*(iz-1)/rz);
        end
    end
    V=repmat(V(:,1),[1,n]);
    V=V*(nx*nys*nzs);
    % V=V/ry/rz;
    m=sparse(I(:),J(:),V(:),n,n);
end