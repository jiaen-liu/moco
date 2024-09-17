function d=pha_crct_epi(d,dblpo,ro_pol,ref_pol,ord,pc_by_echo)
    if isempty(dblpo)
        d=d;
        return;
    end
    if nargin<3
        ro_pol=zeros(necho,1);
        ro_pol(1:2:end)=1;
        ref_pol=1;
    end
    if nargin<5
        ord=2;
    end
    if nargin<6
        pc_by_echo=0;
    end
    [nr,necho,ns,nch,ntr]=size(d);
    if necho<3
        pc_by_echo=0;
    end
    if pc_by_echo
        pha=zeros(nr,necho);
        c=zeros(ord+1,necho);
        for iepc=1:necho
            if iepc==1
                dblpl_tmp=dblpo(:,1:3,:,:);
                [pha(:,iepc),c(:,iepc)]=dpOddEven(dblpl_tmp,ord);
            elseif iepc==necho
                dblpl_tmp=dblpo(:,necho-2:necho,:,:);
                [pha(:,iepc),c(:,iepc)]=dpOddEven(dblpl_tmp,ord);
                if mod(iepc-2,2)==0
                    pha(:,iepc)=-pha(:,iepc);
                    c(:,iepc)=-c(:,iepc);
                end
            else
                dblpl_tmp=dblpo(:,iepc-1:iepc+1,:,:);
                [pha(:,iepc),c(:,iepc)]=dpOddEven(dblpl_tmp,ord);
                if mod(iepc-1,2)==0
                    pha(:,iepc)=-pha(:,iepc);
                    c(:,iepc)=-c(:,iepc);
                end
            end
        end
        c_median=median(c,2);
        x=([0:nr-1].'-floor(nr/2))*1/nr;
        A=x.^[0:ord];
        ref_pha=A*c_median;
    else
        ref_pha=dpOddEven(dblpo,ord);
    end
    d(:,find(ro_pol==ref_pol(1)),:,:,:)=d(:,find(ro_pol==ref_pol(1)),:,:,:)./...
        exp(1i*ref_pha/2);
    d(:,find(ro_pol~=ref_pol(1)),:,:,:)=d(:,find(ro_pol~=ref_pol(1)),:,:,:).*...
        exp(1i*ref_pha/2);

end
