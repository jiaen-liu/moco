function [uniq,ind_uniq,nuniq]=unique_dim(A,dim)
    sia=size(A);
    if numel(sia)<dim
        error('*** The target dimension is beyond the data! ***');
    end
    idx_permute=[1:numel(sia)];
    idx_permute(dim)=1;
    idx_permute(1)=dim;
    
    A=permute(A,idx_permute);
% $$$     B=permute(B,idx_permute);
    si=size(A);
    n=numel(A);
    n1=si(1);
    A=reshape(A,[n1,n/n1]);
% $$$     B=reshape(B,[n1,n/n1]);

    uniq=zeros(n1,n/n1);
    ind_uniq=zeros(n1,1);
    
    nuniq=0;
    nleft=n1;
    mask_left=true(n1,1);
    Aleft=A;
    for i=1:n1
        nleft=nnz(mask_left);
        if nleft==0
            break;
        end
        ind_left=find(mask_left);
        nuniq=nuniq+1;
        uniq(nuniq,:)=Aleft(ind_left(1),:);
        mask_left(ind_left(1))=0;
        ind_uniq(ind_left(1))=i;
        for j=2:nleft
            if isequal(Aleft(ind_left(j),:),Aleft(ind_left(1),:))
                mask_left(ind_left(j))=0;
                ind_uniq(ind_left(j))=i;
            end
        end
        
    end
    if nuniq<n1
        uniq(nuniq+1:n1,:)=[];
    end
    uniq=reshape(uniq,[nuniq,si(2:end)]);
    uniq=permute(uniq,idx_permute);
end