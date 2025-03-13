function y=file_addext(fn,string,n)
% n defines the number of '.' in ext
    if nargin<3
        n=1;
    end
    [fp,name,ext]=fileparts(fn);
    if contains(name,'.') && n>1
        idx=strfind(name,'.');
        if n>(length(idx)+1)
            n=length(idx)+1;
        end
        ext=[name(idx(end-(n-2)):end),ext];
        name=name(1:idx(end-(n-2))-1);
    end
    y=fullfile(fp,[name string ext]);
end