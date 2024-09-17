% MODIFICATION HISTORY:
%    2022-07-13, JAdZ
%        Bugfix: The original version of this tossed all characters in tFree
%        that were not 't' or 'f'.
%    2022-08-16, JAdZ
%        Removed commented-out old code.

function h=siemens_asc_header(file)
    if ~ischar(file)
        file=mid2filename(file);
    end
    prot=eval_twix_hdr(file);
    if iscell(prot)
        prot=prot{end};
    end
    h=prot.MeasYaps;
    if isfield(h,'sWiPMemBlock') && ~isfield(h,'sWipMemBlock')
        h.sWipMemBlock=h.sWiPMemBlock;
        h=rmfield(h,'sWiPMemBlock');
    end
    h.sWipMemBlock.tFree=h.sWipMemBlock.tFree(find(h.sWipMemBlock.tFree ~= '"'));
    nalfree=length(h.sWipMemBlock.alFree);
    alFree=zeros(64,1,'int32');
    for i=1:nalfree
        if ~isempty(h.sWipMemBlock.alFree{i})
            alFree(i)=h.sWipMemBlock.alFree{i};
        end
    end
    h.sWipMemBlock.alFree=alFree;
    
    nadfree=length(h.sWipMemBlock.adFree);
    adFree=zeros(16,1,'double');
    for i=1:nadfree
        if ~isempty(h.sWipMemBlock.adFree{i})
            adFree(i)=h.sWipMemBlock.adFree{i};
        end
    end
    h.sWipMemBlock.adFree=adFree;
end
