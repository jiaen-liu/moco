function y=steRegress(data,r,ord,mask)
    if nargin<3
        ord=floor(r/2-1e-5);
    end
    si=size(data);
    nt=si(end);
    if nargin<4
        mask=true(nt,1);
    end
    n=numel(data)/nt;
    data=reshape(data,[n,nt]);
    x=mtxHarmRegr(nt,1,r,ord);
    y=zeros(n,nt);
    for i=1:n
        if isreal(data)
            [b,~,resi,~,stat]=regress(col(data(i,mask(:))),x(mask(:),:));
            y(i,:)=data(i,:)-(x(:,2:end)*b(2:end))';
        else
            [b,~,resi,~,stat]=regress(real(col(data(i,mask(:)))),x(mask(:),:));
            tmpr=real(data(i,:))-(x(:,2:end)*b(2:end))';
            [b,~,resi,~,stat]=regress(imag(col(data(i,mask(:)))),x(mask(:),:));
            tmpi=imag(data(i,:))-(x(:,2:end)*b(2:end))';
            y(i,:)=tmpr+1i*tmpi;
        end
    end
    y=reshape(y,si);
end