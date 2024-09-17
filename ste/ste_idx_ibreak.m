function idx=ste_idx_ibreak(sorted_ste)
    ishot_break=sorted_ste.ishot_break;
    nparv_cyc=sorted_ste.nparv_cyc;
    if nparv_cyc*sorted_ste.nvf<sorted_ste.nv
        nv=nparv_cyc*(sorted_ste.nvf+1);
    end
    idx=zeros(nv,1);
    idx_shot=sorted_ste.idx_shot(2,:);

    if isstruct(ishot_break)
        fnm=fieldnames(ishot_break);
        n=length(fnm);
        nc=0;
        for i=1:n
            var=getfield(ishot_break,fnm{i});
            if ~isempty(var)
                nc=nc+1;
                nseg=size(var,2);
                for j=1:nseg
                    idx(idx_shot(var(1,j)+1):idx_shot(var(2,j)+1))=nc;
                end
            else
                continue;
            end
        end
    elseif isvector(ishot_break)
        nc=1+length(ishot_break);
        for i=1:nc
            if i==1
                idx(idx_shot(1):idx_shot(ishot_break(1)+1))=i;
            elseif i==nc && sorted_ste.nshot>ishot_break(i)+1
                idx(idx_shot(ishot_break(i)+1):end)=i;
            else
                idx(idx_shot(ishot_break(i-1)+1):idx_shot(ishot_break(i)+1))=i;
            end
        end
    end        
end