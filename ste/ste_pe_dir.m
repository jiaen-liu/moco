function pe_dir=ste_pe_dir(data)
    ncontr=data.pe.ncontr;
    necho=numel(data.te);
    necho_contr=round(necho/ncontr);
    ky=data.pe.k(1,1+necho:2*necho);
    if max(abs(ky))==0
        ky=data.pe.k(1,1+2*necho:3*necho);
    end
    pe_dir=zeros(1,ncontr);
    for i=1:ncontr
        pe_dir(i)=sign(ky(necho_contr*(i-1)+2)-ky(necho_contr*(i-1)+1));
    end
end