function k=genk(para)
    ntr=para.n_main_tr;
    np=para.np;
    ns=para.n_slices;
    npar=para.n_partitions;
    s1=para.sense_rate_p;
    s2=para.sense_rate_s;
    nref=para.n_refs(1);
    necho_intl=para.np/para.n_interleaves/s1;
    necho_shot=necho_intl+nref;
    
    nintl=para.n_interleaves;
    k=zeros(2,necho_shot,ntr);
    for i=1:ntr
        pe1=-floor(np/2)+[0:necho_intl-1]*nintl*s1+...
            mod(i-1,nintl)*s1;
        k(1,1:floor(necho_intl/2)+1,i)=pe1(1:floor(necho_intl/2)+1);
        if nref>0
            k(1,floor(necho_intl/2)+2:floor(necho_intl/2)+1+nref,i)=...
                pe1(floor(necho_intl/2)+1);
        end
        if necho_intl>2
            k(1,floor(necho_intl/2)+1+nref+1:end,i)=...
                pe1(floor(necho_intl/2)+2:end);
        end
        k(2,:,i)=-floor(npar/2)+...
                 mod(floor((i-1)/nintl)*s2,npar);
    end
end