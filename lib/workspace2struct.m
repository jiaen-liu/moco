function y=workspace2struct(s)
% s is passed as who()
    n=length(s);
    y=struct;
    for i=1:n
        y=setfield(y,s{i},evalin('caller',s{i}));
    end
end