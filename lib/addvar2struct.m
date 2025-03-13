function y=addvar2struct(x,varargin)
    nvararg=numel(varargin);
    m=zeros(nvararg,1);
    v=cell(nvararg,1);
    y=x;
    for i=1:nvararg
        if ischar(varargin{i})
            y=setfield(y,varargin{i},evalin('caller',varargin{i}));
        end
    end
end