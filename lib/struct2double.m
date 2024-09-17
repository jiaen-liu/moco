function y=struct2double(x)
    if numel(x)>1
        for i=1:numel(x)
            x(i)=struct2double(x(i));
        end
        y=x;
        return;
    end
    names=fieldnames(x);
    n=numel(names);
    for i=1:n
        if isnumeric(x.(names{i})) && ~isa(x.(names{i}),'double')
            x.(names{i})=double(x.(names{i}));
        elseif isstruct(x.(names{i}))
            x.(names{i})=struct2double(x.(names{i}));
        end
    end
    y=x;
end