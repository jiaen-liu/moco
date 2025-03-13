function f=tukey_fraction(si,r)
    nd=length(si);
    ws=zeros(nd,1);
    for i=1:nd
        ws(i)=sum(tukeywin(si(i),r).^2);
    end
    f=prod(ws)/prod(si);
end