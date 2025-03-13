function [y,mdh]=readmdh(mid,acqType)
    if nargin<2
        acqType='all';
    end
    para=extract_para(mid);
    if ~ischar(mid)
        fnmdh=get_file_filter('.',['MID*',num2str(mid),'.mdh']);
    else
        error('*** MID is needed! ***');
    end
    
    nch=para.n_channels;
    % dummy shots were already discarded in sort_siemens.
    ntr=para.n_tr-para.n_dummy_tr;
    ns=para.n_slices;
    nnav=0;
    if para.b_nav_en
        nnav = conditional(para.nav_type==0,1,para.nav_type);
    end
    necho_main=0;
    for icontr=1:para.n_contrasts
        necho_main=necho_main+...
            (para.np/para.n_interleaves/para.sense_rate_p+para.n_refs(icontr));
    end
    % total echos in one shot
    nechot = nnav+para.n_echo_steref+necho_main;
    % read mdh;
    if contains(para.idea_v,'VB','IgnoreCase',true)
        mdh=defMDH17();
        nByteMdh=size_of(mdh);
        size_y=[nByteMdh,nch,nechot,ns,ntr];
        y=read_raw(fnmdh,size_y,'uint8',1);
    elseif contains(para.idea_v,'VD','IgnoreCase',true) || ...
            contains(para.idea_v,'VE','IgnoreCase',true)
        mdh=defMDH11();
        nByteMdh=size_of(mdh);
        size_y=[nByteMdh,1,nechot,ns,ntr];
        y=read_raw(fnmdh,size_y,'uint8',1);
    end
    % for cmrr 7T data (ve12u), it seems slice_order doesn't lead to 
    % interleaved slice order
    if para.dimen==2 && para.slice_order==4 && ~strcmp(acqType,'noise') && ...
                  ~contains(para.idea_v,'VE12U','IgnoreCase',true)
        idx_slice=siemens_slice_order(ns);
        y=y(:,:,:,idx_slice,:);
    end
    y=reshape(y,[size_y(1:3),ns*ntr]);
    if ~strcmp(acqType,'all')
        if strcmp(acqType,'noise')
            y=y(:,:,:,1);
        end
    end
end
