function y=readSortSiem(mid,varargin)
% $$$ % example:
% $$$     d=readSortSiem(mid,'main',1);
% $$$     d=readSortSiem(mid,'noise',1);
% $$$     % read all
% $$$     d=readSortSiem(mid);
% $$$     % read TR by idxtr
% $$$     d=readSortSiem(mid,'idxtr',[2,5]);
% Modifiction history:
%    2022-10-04 JAdZ
%        Now calling sort_siemens with the -filemode 0666 option to allow all to
%        read and overwrite.
% parse input
    p=inputParser;
    keyNoise=0;
    keyBlipoff=0;
    keyMain=0;
    idxTR=[];
    keyNav=0;
    addParameter(p,'noise',keyNoise,@isnumeric);
    addParameter(p,'blipoff',keyBlipoff,@isnumeric);
    addParameter(p,'main',keyMain,@isnumeric);
    addParameter(p,'nav',keyNav,@isnumeric);
    addParameter(p,'idxtr',idxTR,@isnumeric);
    p.parse(varargin{:});
    keyNoise=p.Results.noise;
    keyBlipoff=p.Results.blipoff;
    keyMain=p.Results.main;
    keyNav=p.Results.nav;
    idxTR=p.Results.idxtr;
    para=extract_para(mid);
    % read files with different names
    if keyNav
        ff=['MID*',num2str(mid),...
            '.nav.svd'];
    else
        ff=['MID*',num2str(mid),...
            '.raw.svd'];
    end
    fn=get_file_filter('.',ff);
    if isempty(fn)
        cmd=['sort_siemens -filemode 0666 ' num2str(mid)];
        if system(cmd)~=0
            error(['*** ',cmd,' was not successful! ***']);
        end
    end
    if isempty(idxTR)
        % idxTR is determined based on acqusition type
        
        % determine the number of volume TRs for 
        % each acqusition type
        % Dummy shots are removed by sort_siemens
        % 'idxTRDummy','n_dummy_tr' is not needed after noise
        strAcq={'idxTRNoise','n_noise_shots';...
                'idxTRBlipoff','n_blipoff_tr';...
                'idxTRnCycVte','n_nvarte_tr';...
                'idxTRMain','n_main_tr'};
        i1=1;
        i2=0;
        nAcq=size(strAcq,1);
        for i=1:nAcq
            i2=i2+para.(strAcq{i,2});
            if i2<i1
                eval([strAcq{i,1},...
                     '=[];']);
            else
                
                eval([strAcq{i,1},...
                     '=[i1,i2];']);
            end
            i1=i1+para.(strAcq{i,2});
        end
        if keyNoise
            y=read_data(fn,'partial',[5,idxTRNoise]);
            % for noise scan
            % data in the third dimention 
            % contains zeros
            nn=size(y,3);
            idx=zeros(nn,1);
            for ii=1:nn
                if min(col(abs(y(:,:,ii,:,:))))>0
                    idx(ii)=1;
                end
            end
            y=y(:,:,find(idx),:,:);
        elseif keyBlipoff
            y=read_data(fn,'partial',[5,idxTRBlipoff]);
        elseif keyMain
            y=read_data(fn,'partial',[5,idxTRMain]);
        else
            % if no key words are provided
            % read all data
            y=read_data(fn);
        end
    elseif numel(idxTR)==1
        y=read_data(fn,'partial',[5,idxTR(1),idxTR(1)]);
    elseif numel(idxTR)==2
        y=read_data(fn,'partial',[5,idxTR(1),idxTR(2)]);
    else
        error('*** Error about the TR index! ***');
    end
    %% Multislicemode is actually not used by amri_epi
    % instead, use ucmode or slice_order to find the slice ordering
    if para.slice_order ~= 1 && para.slice_order ~= 4
        error('*** Slice ordering is not supported! ***');
    end
    % for cmrr 7T data (ve12u), it seems slice_order doesn't lead to 
    % interleaved slice order
    if para.dimen==2 && para.slice_order==4 && size(y,3)>1 && ...
                  ~contains(para.idea_v,'VE12U','IgnoreCase',true)
        idx_slice=siemens_slice_order(para.n_slices);
        y=y(:,:,idx_slice,:,:);
    end
end
