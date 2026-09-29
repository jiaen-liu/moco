function par=readFilePar(fn)
% $$$     file example
% $$$     // comment
% $$$     $a 1
% $$$     $b 2
% $$$     $ninv 3 4 // comment
% $$$     $c 4
% $$$     $d %Name
% $$$     $e #[0.1:0.1:0.5].'*1e-3
% $$$     gives
    
% $$$     par = 
% $$$ 
% $$$   struct with fields:
% $$$ 
% $$$        a: 1
% $$$        b: 2
% $$$     ninv: [3 4]
% $$$        c: 4
% $$$        d: 'Name'
% $$$        e: [0.1:0.1:0.5].'*1e-3

    %% important
    % space is counted as a separator. Avoid using space if not necessary
    if ~exist(fn)
        error('*** File doesn''t exist! ***');
    end
    fid=fopen(fn);
    c=textscan(fid,'%s','CommentStyle','//');
    c=c{1};
    n=length(c);
    par=struct;
    for i=1:n
        if strcmp(c{i}(1),'$')
            varName=c{i}(2:end);
            par=setfield(par,varName,[]);
        elseif strcmp(c{i}(1),'%')
            par.(varName)=c{i}(2:end);
        elseif strcmp(c{i}(1),'#')
            par.(varName)=eval(c{i}(2:end));
        else
            % number
            par.(varName)=[par.(varName),str2num(c{i})];
        end
    end
    fclose(fid);
end
