function [arr]=sense_ref_philips(mid,varargin)
    type=1;
    apodiz=0.25;
    corr_half_fov=0;
    p=inputParser;
    dimen=3;
    addParameter(p,'type',type,@isnumeric);
    addParameter(p,'corr_half_fov',corr_half_fov,@isnumeric);
    addParameter(p,'dimen',dimen,@isnumeric);
    p.parse(varargin{:});
    type=p.Results.type;
    corr_half_fov=p.Results.corr_half_fov;
    dimen=p.Results.dimen;
    r=MRecon(mid2filename_philips(mid));
    % r.Parameter.Parameter2Read.typ=type;
    r.ReadData;
    % preprocessing
    r.DcOffsetCorrection;
    r.PDACorrection;
    r.RandomPhaseCorrection;
    r.MeasPhaseCorrection;
    label=r.Parameter.Labels.Index;
    dataall=r.Data{1,1};

    idx_typ1=find(label.typ==type);
    loca=label.loca(idx_typ1);
    chanall=label.chan(idx_typ1);
    kyall=label.ky(idx_typ1);
    kzall=label.kz(idx_typ1);
    % idx=find(loca==0&chan>1);
    if dimen==2
        ns=length(unique(loca));
    else
        ns=1;
    end
    arr=cell(ns,1);
    for ind_sl=1:ns
        idx=find(loca==ind_sl-1);
        ky=kyall(idx);
        kz=kzall(idx);
        chan=chanall(idx);
        data=dataall(:,idx);

        nch=length(unique(chan));
        
        ky=ky(1:nch:end);
        kz=kz(1:nch:end);
        
        [nr,ntmp]=size(data);
        nkr=nr;
        nk=ntmp/nch;
        data=reshape(data,[nr,nch,nk]);
        data=permute(data,[1,3,2]);
        
        k=[ky,kz].';

        % correct half fov shift
        idx_rev=mod(k(1,:),2)==1;
        data(:,idx_rev,:)=-data(:,idx_rev,:);
        idx_rev=mod(k(2,:),2)==1;
        data(:,idx_rev,:)=-data(:,idx_rev,:);
        
        nkp=r.Parameter.Encoding.KyRange(1,2)-r.Parameter.Encoding.KyRange(1,1)+1;
        if isempty(r.Parameter.Encoding.KzRange)
            nks=1;
        else
            nks=r.Parameter.Encoding.KzRange(1,2)-r.Parameter.Encoding.KzRange(1,1)+1;
        end
        arr{ind_sl}=zeros(nkr,nkp,nks,nch);
        
        for ik=1:nk
            ip=k(1,ik)+floor(nkp/2)+1;
            is=k(2,ik)+floor(nks/2)+1;
            arr{ind_sl}(:,ip,is,:)=data(:,ik,:);
        end
        % correct readout oversampling
        arr{ind_sl}=fftmr(arr{ind_sl},-1,1);
        
        arr{ind_sl}=arr{ind_sl}(idx_truncate(nr,nr/2),:,:,:);
        
        arr{ind_sl}=fftmr(arr{ind_sl},1,1);
        % reduce gibbs ring effect
        if dimen==3
            arr{ind_sl}=arr{ind_sl}.*win_tukey(nr/2,apodiz).*...
                reshape(win_tukey(nkp,apodiz),[1,nkp]).*...
                reshape(win_tukey(nks,apodiz),[1,1,nks]);
        else
            arr{ind_sl}=arr{ind_sl}.*win_tukey(nr/2,apodiz).*...
                reshape(win_tukey(nkp,apodiz),[1,nkp]);
        end
        arr{ind_sl}=fftmr(arr{ind_sl},-1,[1,2,3]);
    end
    if ns>1
        arr_tmp=zeros([size(arr{1}),ns]);
        for ind_sl=1:ns
            arr_tmp(:,:,:,:,ind_sl)=arr{ind_sl};
        end
        arr_tmp=permute(arr_tmp,[1,2,5,4,3]);
        arr=arr_tmp;
    end
    
end