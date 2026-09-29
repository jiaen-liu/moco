function y=combine_dim(d,dim)
    si=size(d);
    ndim=length(si);
    ndim_c=length(dim);
    
    
    if dim(1)==1
        si_c(1)=prod(si(dim));
        if dim(end)<ndim
            si_c=[prod(si(dim)),col(si(dim(end)+1:end)).'];
        elseif dim(end)==ndim
            si_c=[prod(si(dim)),1];
        else
            error('*** The requested dimension is more than the data dimension! ***');
        end
    else
        si_c=zeros(1,ndim-ndim_c+1);
        si_c(1:dim(1)-1)=si(1:dim(1)-1);
        si_c(dim(1))=prod(si(dim));
        if dim(end)~=ndim
            si_c(dim(1)+1:end)=si(dim(end)+1:end);
        end
    end
    y=reshape(d,si_c);
end
