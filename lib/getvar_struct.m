function getvar_struct(astruct,varargin)
    nvarargin=numel(varargin);
    if nvarargin > 0
        fldn=varargin;
    else
        fldn=fieldnames(astruct);
    end
    nfldn=numel(fldn);
    for i=1:nfldn
        if isfield(astruct,fldn{i})
            assignin('caller',fldn{i},getfield(astruct,fldn{i}));
        else
            assignin('caller',fldn{i},[]);
        end
    end
end