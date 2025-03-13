function [data,label]=read_raw_philips(mid,varargin)
% use reconframe
% type: 1: image data, 3, phase correction, 5, noise data
% mix: 0: regular data, 1: 3d navigator
    type=1;
    mix=0;
    corr_half_fov=0;
    p=inputParser;
    addParameter(p,'type',type,@isnumeric);
    addParameter(p,'mix',mix,@isnumeric);
    addParameter(p,'corr_half_fov',corr_half_fov,@isnumeric);
    p.parse(varargin{:});
    type=p.Results.type;
    mix=p.Results.mix;
    corr_half_fov=p.Results.corr_half_fov;
    r=MRecon(mid2filename_philips(mid));
    r.Parameter.Parameter2Read.typ=type;
    r.Parameter.Parameter2Read.mix=mix;
    r.ReadData;
    if mix==0
        r.DcOffsetCorrection;
    end
    r.PDACorrection;
    if mix~=1
        r.RandomPhaseCorrection;
    end
    r.MeasPhaseCorrection;
    data=r.Data;
    % label=r.DataProperties.labels.labels;
    label=r.Parameter.Labels.Index;
    para=[];
    try
        para=extract_para_philips(mid);
    catch me
        if strcmp(me.message,'The parameter does not exist')
            warning(['*** Parameter not saved in the raw data. Default ones will be used and can be wrong! ***']);
            para.isgre=1;
        end
    end
    nch=size(r.Parameter.Labels.CoilNrs,1);
    nx=size(data,1);
    % para.n_channels;
    names=fieldnames(label);
    nlabel=length(names);
    idx_label=find(label.typ(:)==type & label.mix(:)==mix);
    for i=1:nlabel
        label.(names{i})=col(label.(names{i})(idx_label));
    end
    switch type
      case 1
        if isempty(r.Data)
            data=[];
            label=[];
            return;
        end
        if corr_half_fov
            if para.isgre
                idx_rev=mod(label.ky(:),2)==1;
                data(:,idx_rev)=-data(:,idx_rev);
            end
            idx_rev=mod(label.kz(:),2)==1;
            data(:,idx_rev)=-data(:,idx_rev);
        end
      case 3
      case 5
      case 7
        % phase navigator
        nnav=para.nav1d_enable;
        if nnav>1
            data=reshape(data,[nx,nch,nnav,numel(data)/nx/nch/nnav]);
            % in current reconframe, the direction of the second phase
            % navigator is always flipped. The code below is to undo
            % the first flip when the gradients have the same polarity
            if sign(r.Parameter.GetValue('GR::phnav_0_->str'))==...
                    sign(r.Parameter.GetValue('GR::phnav_1_->str'))
                data(:,:,2,:)=flipdim(data(:,:,2,:),1);
            end
        end
      otherwise
    end
    for i=1:nlabel
        tmp=label.(names{i});
        label.(names{i})=double(tmp(1:nch:end));
    end
    % With 34 channels, the last two channels are tossed out for EPI phase correction specific to the 7T
    if abs(para.frequency/42.58e6-7)<0.5
        if nch>32
            if nch~=34
                error('*** There should be 32 or 34 channels! ***');
            end
            s1=size(data,1);
            data=reshape(data,[s1,nch,numel(data)/...
                               s1/nch]);
            channel_id=r.Parameter.Labels.CoilNrs(:,1);
            data=data(:,channel_id>1,:);
            data=reshape(data,[s1,numel(data)/s1]);
        end
    end
end