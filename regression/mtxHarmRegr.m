function X=mtxHarmRegr(n,dt,r,ord,mask)
    if nargin<5
        mask=true(n,1);
    end
    t=[0:n-1].'*dt;
    if ord~=r/2
        X=zeros(n,1+ord*2);
        X(:,1)=1;
        for i=1:ord
            X(:,2*i)=cos(2*pi/r*t*i);
            X(:,2*i+1)=sin(2*pi/r*t*i);
        end
    else
        X=zeros(n,ord*2);
        X(:,1)=1;
        for i=1:ord-1
            X(:,2*i)=cos(2*pi/r*t*i);
            X(:,2*i+1)=sin(2*pi/r*t*i);
        end
        X(:,2*ord)=cos(2*pi/r*t*ord);
    end
    X=X(find(mask),:);
end