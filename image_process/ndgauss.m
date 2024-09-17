function h=ndgauss(siz,sig)
    nd=numel(siz);
    if numel(sig)~=nd && numel(sig)==1
        sig=ones(nd,1)*sig;
    elseif numel(sig)~=nd
        error('The dimension of sigma should be equal to that of size');
    end
    
    siztmp=zeros(2,length(siz));
    siztmp(1,:)=-floor(siz/2);
    siztmp(2,:)=floor((siz-1)/2);
    siz=siztmp;
    x=cell(nd,1);
    
    expression='[';
    for i=1:nd
        if i<nd
            expression=[expression 'x{' int2str(i) '},'];
        else
            expression=[expression 'x{' int2str(i) '}]'];
        end
    end
    expression=[expression '=ndgrid('];
    for i=1:nd
        if i<nd
            expression=[expression 'siz(1,' int2str(i) '):siz(2,' int2str(i) ...
                        '),'];

        else
            expression=[expression 'siz(1,' int2str(i) '):siz(2,' int2str(i) ...
                        '));'];
        end
    end
    eval(expression);
    A=x{1}.*x{1}/2/sig(1)^2;
    if nd>1
        for i=2:nd
            A=A+x{i}.*x{i}/2/sig(i)^2;
        end
    end
    h=exp(-A);
    h=h/sum(h(:));
end