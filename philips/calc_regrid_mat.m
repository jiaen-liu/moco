function A=calc_regrid_mat(nus,nr_os,apodiz)
    nk=length(nus);
    np=min(nk,nr_os);
    dnus=max(abs(nus(2:end)-nus(1:end-1)));
    np=floor(nr_os/dnus);
    c=100;
    while true
        pos=[0:np-1]-floor(np/2);
        % nus in unit of 1/sampled_FOV
        k=nus(:)*(2*pi/nr_os);
        A=exp(1i*k.*pos);
        c=cond(A);
        % 50 is emperical
        if c>50
            np=np-1;
        else
            break;
        end
    end
    A=inv(A'*A)*A';
    if apodiz ~= -1
        % apodization
        k_min = min(k(:));
        k_max = max(k(:));
        w = tukeyWinNUGrid(nk, apodiz, (k(:)+pi)/(2*pi));
        A=A*diag(w);
    end
end