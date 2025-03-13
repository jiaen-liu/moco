% MODIFICATION HISTORY:
%    2022-07-11, JAdZ
%        The wipmem_l array is now expanded with alTD content.
%    2022-07-12, JAdZ
%        Support for dense tFree use near the end of the string.
%    2022-07-13, JAdZ
%        Changes to this header, no functional changes.
%    2022-08-16, Jiaen Liu
%        Add condition for alTD, which doesn't exit in some idea versions.

function wip=amri_epi_wipmem(header)
    if ~isstruct(header)
        header=siemens_asc_header(header);
    end
    wipdef=define_amri_epi_wipmem();
    wipb_tmp=wipbool_extract(header.sWipMemBlock.tFree);
    wipmem_b=cast2struct(wipb_tmp,wipdef.wipmem_b);
    wipmem_d=cast2struct(header.sWipMemBlock.adFree,wipdef.wipmem_d);
	% 2022-07-11 JAdZ - make sure that we have 256 bytes in wipl, then append alTD(33:..)
	wipl=typecast(header.sWipMemBlock.alFree,'uint8');
	while size_of(wipl) < 256
		wipl(end+1)=uint8(0);
	end
    if isfield(header,'alTD')
        for i=33:numel(header.alTD)
            if ~isempty(header.alTD{i})
                wipl(end+1:end+4)=typecast(uint32(cell2mat(header.alTD(i))),'uint8');
            else 
                wipl(end+1:end+4)=typecast(uint32(0),'uint8');
            end
        end
    end
	while size_of(wipl) < size_of(wipdef.wipmem_l)
		wipl(end+1)=uint8(0);
	end
    wipmem_l=cast2struct(wipl,wipdef.wipmem_l);
    wip=struct;
    wip=pass_var_struct(wip,wipmem_l);
    wip=pass_var_struct(wip,wipmem_d);
    wip=pass_var_struct(wip,wipmem_b);
end
