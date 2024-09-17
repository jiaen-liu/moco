function pe=get_pe_philips(mid,varargin)
% use reconframe
% type: 1: image data, 3, phase correction, 5, noise data
% mix: 0: regular data, 1: 3d navigator
    type=1;
    mix=0;
    p=inputParser;
    addParameter(p,'type',type,@isnumeric);
    addParameter(p,'mix',mix,@isnumeric);
    p.parse(varargin{:});
    type=p.Results.type;
    mix=p.Results.mix;
    r=MRecon(mid2filename_philips(mid));
    r.Parameter.Parameter2Read.typ=type;
    r.Parameter.Parameter2Read.mix=mix;
    label=r.Parameter.Labels.Index;
    para=extract_para_philips(mid);
    nch=para.n_channels;
    ky=label.ky(label.typ==type&label.mix==mix);
    kz=label.kz(label.typ==type&label.mix==mix);
    ky=ky(1:nch:end);
    kz=kz(1:nch:end);
    s1=para.sense_rate_p;
    s2=para.sense_rate_s;
    pe=[ky(:).'*s1;kz(:).'*s2];
    pe=reshape(pe,[2,para.nte_contr*para.epi_factor,...
                   numel(pe)/2/para.nte_contr/para.epi_factor]);
end