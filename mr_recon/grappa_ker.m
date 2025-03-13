function [kp,ind_uniq]=grappa_ker(ry,rz,dz,nky,nkz)
% this function calculate the neighborhood source points for 2d-grappa-caipi reconstruction
% ry rz are acceleration ratio in y and z direction
% dz is shift of kyz line for every kz
% nky and nkz are the grappa kernel size


% $$$     example of ryxrz=3x3, dz=0
% $$$     o x x o * * o * * o 
% $$$     x x x * * * * * * *
% $$$     x x x * * * * * * *
% $$$     o * * o * * o * * o
% $$$     example of ryxrz=3x3, dz=2
% $$$     o x x o * * o * * o 
% $$$     x x x * * * * * * *
% $$$     x x x * * * * * * *
% $$$     * * o * * o * * o *
% determine the number of positions: nr (x marker)
% determine the neighbor source points for x 
    nr=ry*rz;
    kp=cell(nr-1,1);
    for i=1:nr-1
        [iy,iz]=ind2sub([ry,rz],i+1);
        % location in k space
        kp{i}.iy=iy;
        kp{i}.iz=iz;
        % ky location of source points
        if iy-1<=ry+1-iy
            % closer to ky=1
            % 1, 1-ry,1-ry*2,1-ry*3...
            % ry+1, ry*2+1,ry*3+1...
            ys_odd=[1:-ry:1-(ceil(nky/2)-1)*ry];
            ys_even=[1+ry:ry:1+floor(nky/2)*ry];
        else
            % closer to ky=ry+1
            % ry+1, ry*2+1,ry*3+1...            
            % 1, 1-ry,1-ry*2,1-ry*3...
            ys_odd=[1+ry:ry:1+ceil(nky/2)*ry];
            ys_even=[1:-ry:1-(floor(nky/2)-1)*ry];
        end
        ys=comb_interleave(ys_odd,ys_even);
        % define the first grid
        if iz-(dz/ry)*(iy-1)-1<=0
            % iy>=(ry/dz)*(iz-1)+1
            z1=1-rz;
            z2=1;
        else
            z1=1;
            z2=1+rz;
        end
        if iz-((z1+z2)/2+(dz/ry)*(iy-1))<=0
            % closer to the left edge
            zs_odd=[z1:-rz:z1-(ceil(nkz/2)-1)*rz];
            zs_even=[z1+rz:rz:z1+floor(nkz/2)*rz];
        else
            % closer to the right edge
            zs_odd=[z1+rz:rz:z1+ceil(nkz/2)*rz];
            zs_even=[z1:-rz:z1-(floor(nkz/2)-1)*rz];
        end
        zs=comb_interleave(zs_odd,zs_even);
        % convert into grid
        [ys,zs]=ndgrid(ys,zs);
        ys=ys(:);
        zs=zs(:);
        % **************** %
        % normalized position of the kernel in unshifted coordinate
        kp{i}.ysn=(ys-1)/ry;
        kp{i}.zsn=(zs-1)/rz;
        % **************** %
        % real position relative to each position
        % account for the kz shift
        kp{i}.ys=ys-iy;
        zs=(ys-1)/ry*dz+zs; 
        kp{i}.zs=zs-iz;
        
        
    end
    % find unique kernels
    yzsn=zeros(2*(nky*nkz),nr-1);
    for i=1:nr-1
        yzsn(1:nky*nkz,i)=kp{i}.ysn;
        yzsn(nky*nkz+1:end,i)=kp{i}.zsn;
    end
    [~,ind_uniq,nuniq]=unique_dim(yzsn,2);
    for i=1:nr-1
        kp{i}.ind_uniq=ind_uniq(i);
    end
end



