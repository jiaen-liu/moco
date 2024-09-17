function save_mat(fn,varargin)
% use instead of matlab's save: no overwrite, check file size>2GB
    overwrite=0;
    read4all=0;
    p=inputParser;
    nvar=length(varargin);
    bstrcmp=strcmp(varargin,'overwrite');
    if any(bstrcmp)
        idx=find(bstrcmp);
        addParameter(p,'overwrite',overwrite,@isnumeric);
        p.parse(varargin{idx:idx+1});
        overwrite=p.Results.overwrite;
        nvar=nvar-2;
    end
    bstrcmp=strcmp(varargin,'read4all');
    if any(bstrcmp)
        idx=find(bstrcmp);
        addParameter(p,'read4all',read4all,@isnumeric);
        p.parse(varargin{idx:idx+1});
        read4all=p.Results.read4all;
        nvar=nvar-2;
    end    
    if isfile(fn) && ~overwrite
        disp(['*** File ''' fn ''' exist! ***']);
        fn_test=fn;
        for i=1:100
            if exist(fn_test)==2
                fn_test=file_addext(fn,['_' num2str(i)]);
            else
                break;
            end
        end
        fn=fn_test;
        disp(['*** It will be saved as ''' fn ''' ***']);
    end
    
% $$$     if overwrite==1 || strcmp(varargin{nvar-1},'overwrite')
% $$$         nvar=nvar-2;
% $$$     end
    nbyte=0;
    cmd=['save(''' fn ''''];
    for i=1:nvar
        cmd_info=['whos(''' varargin{i} ''')'];
        info=evalin('caller',cmd_info);
        % tmp=evalin('caller',varargin{i});
        % assignhere(varargin{i},tmp);
        % info=whos(varargin{i});
        nbyte=nbyte+info.bytes;
        cmd=[cmd ',''' varargin{i} ''''];
    end


    if nbyte>2*1024^3*0.95
        % use -7.3
        disp('*** Save in v7.3 ***');
        cmd=[cmd ',''-v7.3'');'];
    else
        cmd=[cmd ');'];
    end
    evalin('caller',cmd);
    if read4all
        file_permission(fn,'+rw','ugo');
    end
end
