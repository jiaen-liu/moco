function y=ste_interp(x,reso,resn)
    si=size(x);
    nxo=si(1);
    nyo=si(2);
    nzo=si(3);

    fovx=reso(1)*nxo;
    fovy=reso(2)*nyo;
    fovz=reso(3)*nzo;
    
    xo=([1:nxo]-(1+nxo)/2)*reso(1);
    yo=([1:nyo]-(1+nyo)/2)*reso(2);
    zo=([1:nzo]-(1+nzo)/2)*reso(3);

    [xo,yo,zo]=ndgrid(xo,yo,zo);

    nxn=round_even(fovx/resn(1));
    nyn=round_even(fovy/resn(2));
    nzn=round_even(fovz/resn(3));

    xn=([1:nxn]-(1+nxn)/2)*resn(1);
    yn=([1:nyn]-(1+nyn)/2)*resn(2);
    zn=([1:nzn]-(1+nzn)/2)*resn(3);

    [xn,yn,zn]=ndgrid(xn,yn,zn);

    coordo=cat(4,xo,yo,zo);
    coordn=cat(4,xn,yn,zn);
    y=interp3_nmat(coordo,coordn,x,'spline');
end