function y=var2struct(varargin)
    narg=numel(varargin);
    m=zeros(narg,1);
    v=cell(narg,1);
    y=struct;
    for i=1:narg
        if ischar(varargin{i})
            y=setfield(y,varargin{i},evalin('caller',varargin{i}));
        end
    end
end