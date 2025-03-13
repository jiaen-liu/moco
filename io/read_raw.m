function y=read_raw(fn,dim,type,part)
    if exist(fn,'file') ~= 2
        error('The file does not exist!');
    end
    fstat=dir(fn);
    if nargin<4
        part=0;
    end
    if part
        nbyte=prod(dim)*size_of(zeros(1,1,type));
    else
        nbyte=fstat.bytes;
    end
    f=fopen(fn,'r');
    ytmp=fread(f,nbyte,'uint8=>uint8');
    fclose(f);
    y=typecast(ytmp,type);
    if ~isempty(dim) 
        y=reshape(y,dim);
    end
end