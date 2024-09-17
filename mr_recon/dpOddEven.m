% 2023-04-28, Jiaen Liu: avoid phase wrapping issue at the expected
% reference point (intercept) 
% 2022-10-26, Jiaen Liu: In case 20221013_2 the navigator showed lots of artifact It was due to the mask is too small, modified [y,c]=polypha1d(dpc,maskDel,ord,maskDel); to 
% [y,c]=polypha1d(dpc,mask,ord,maskDel); and added masking criteria
function [y,c,dpc,dpc_2e]=dpOddEven(d,ord,frac,nofit)
% data should have [nx,necho,nch,nacq] format
% here nacq can be number of slices or 
% any dimension that will be averaged
% data should be in image domain
    dpc_2e=[];
    if nargin<2
        ord=2;
    end
    if nargin<3
        frac=0.3;
    end
    if nargin<4
        nofit=0;
    end
    [nr,necho,nch,nacq]=size(d);
    if nacq>1
        d=combine_dim(d,[3,4]);
    end
    rms=squeeze(mean(abs(d(:,1,:)).^2,3).^0.5);
    mask=mask1d(rms,frac);
    if mean(mask)>0.7
        [rmssort]=sort(rms(:));
        thrd=rmssort(floor(length(rmssort)*0.3));
        mask=rms>thrd;
    end
    idx_mask=find(mask);
    span_mask=idx_mask(end)-idx_mask(1);
    rel=0.5;
    while true
        maskDel=mask1d(rms,frac,rel);
        idx_mask_del=find(maskDel);
        span_mask_del=idx_mask_del(end)-idx_mask_del(1);
        if span_mask_del/span_mask>0.4
            break;
        else
            rel=rel-0.05;
        end
    end
    if necho>2
        dprc=sum((d(:,1:end-2,:).*d(:,3:end,:)).*...
                 exp(-1i*angle(d(:,2:end-1,:))).^2,3);
        dprc(:,2:2:end,:)=conj(dprc(:,2:2:end,:));
        dprc_2e=sum(d(:,1:2:floor(necho/2)*2,:).*...
                    conj(d(:,2:2:floor(necho/2)*2,:)),3);
    elseif necho==2
        dprc=sum(d(:,1,:,:).*...
                 conj(d(:,2,:,:)),3);
    end
    dprc=sum(dprc,2);
    dpc=angle(dprc);
    if necho>2
        dprc_2e=sum(dprc_2e,2);
        dpc_2e=angle(dprc_2e);
    end
    % fit
    if ~nofit
        if necho>2
            % avoid phase wrapping issue at the expected
            % reference point (intercept)
            idx_mask=find(mask);
            n_center=10;
            idx_mask=idx_mask(1:floor(length(idx_mask)/n_center):end);
            n_center=length(idx_mask);
            dif=zeros(n_center,1);
            y=zeros(nr,n_center);
            c=zeros(ord+1,n_center);
            for i=1:n_center
                [y3e,c3e]=polypha1d(dpc,mask,ord,maskDel,...
                                    idx_mask(i));
                y3e=y3e/2;
                c3e=c3e/2;
                [y2e,c2e]=polypha1d(dpc_2e,mask,ord,maskDel,...
                                    idx_mask(i));
                y(:,i)=y3e;
                c(:,i)=c3e;
                dif(i)=mean(abs(y3e(mask)-y2e(mask)));
            end
            [~,imin]=min(abs(dif));
            % if none center point gives reasonable result
            if abs(dif(imin))>pi/2
                dif_sort=sort(dif);
                dif_median=dif_sort(floor(n_center/2));
                idx_median=find(dif==dif_median,1);
                y=y(:,idx_median)-round(dif_median/pi)*pi;
                c=c(:,idx_median);
                c(1)=c(1)-round(dif_median/pi)*pi;
            else
                y=y(:,imin);
                c=c(:,imin);
            end
        else
            [y,c]=polypha1d(dpc,mask,ord,maskDel);
        end
    else
        y=[];
        c=[];
    end
end
