function [idx,mask_main]=idx_channel_philips(para_ref,para,exclude_channel)
    if nargin<3
        exclude_channel=[];
    end
    idx=zeros(length(para.channel_id),1);
    mask_main=zeros(length(para.channel_id),1);
    for i=1:length(para.channel_id)
        tmp=find(para_ref.channel_id==...
                 para.channel_id(i));
        if isempty(tmp) || ismember(para.channel_id(i),exclude_channel)
            idx(i)=-1;
            mask_main(i)=0;
        else
            idx(i)=tmp;
            mask_main(i)=1;
        end
    end
    idx=idx(idx~=-1);
end
