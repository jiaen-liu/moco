function y=cast2struct(d,structDef)
    d=d(:);
    if ~isa(d,'uint8')
        d=typecast(d,'uint8');
    end
    if size_of(d)<size_of(structDef)
        % to be compatible with previous data
        warning('*** The memory sizes do not match! ***');
    end
    nbytedata=size_of(d);
    tags=fieldnames(structDef);
    ntags=length(tags);
    y=structDef;
    ibyte=0;
    for i=1:ntags
        t=class(y.(tags{i}));
        m=size_of(y.(tags{i}));
        if isnumeric(y.(tags{i}))
            y.(tags{i})=typecast(d(ibyte+1:ibyte+m),t);
        elseif isstruct(y.(tags{i}))
            y.(tags{i})=cast2struct(d(ibyte+1:ibyte+m),y.(tags{i}));
        else 
            error('*** Only numerical type is supported in this function! ***');
        end
        ibyte=ibyte+m;
        % to be compatible with previous data
        if ibyte>=nbytedata
            break;
        end
    end
end
