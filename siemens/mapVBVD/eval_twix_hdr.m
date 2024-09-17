function prot = eval_twix_hdr(filename)
% Author: Philipp Ehses (philipp.ehses@tuebingen.mpg.de), Feb/05/2015
% Updated based on the 2018 version mapVBVD by Jiaen Liu, Nov/19/2020
    if ~exist('filename','var') || isempty(filename)
        info = 'Please select binary file to read';
        [fname,pathname]=uigetfile('*.dat',info); 
        if isempty(pathname)
            return
        end
        filename=[pathname fname];
    else
        if ischar(filename) 
            % assume that complete path is given
            if  ~strcmpi(filename(end-3:end),'.dat');
                filename=[filename '.dat'];   %% adds filetype ending to file
            end
        else
            % filename not a string, so assume that it is the MeasID
            measID   = filename;
            filelist = dir('*.dat');
            filesfound = 0;
            for k=1:numel(filelist)
                if regexp(filelist(k).name,['^meas_MID0*' num2str(measID) '.*\.dat'])==1
                    if filesfound == 0
                        filename = filelist(k).name;
                    end
                    filesfound = filesfound+1;
                end
            end
            if filesfound == 0
                error(['File with meas. id ' num2str(measID) ' not found.']);
            elseif filesfound > 1
                disp(['Multiple files with meas. id ' num2str(measID) ' found. Choosing first occurence.']);
            end
        end
    end
    
    fid = fopen(filename,'r','l','US-ASCII');
    % get file size
    fseek(fid,0,'eof');
    fileSize = ftell(fid);
    % start of actual measurement data (sans header)
    fseek(fid,0,'bof');
    
    firstInt  = fread(fid,1,'uint32');
    secondInt = fread(fid,1,'uint32');

    % lazy software version check (VB or VD?)
    if and(firstInt < 10000, secondInt <= 64)
        version = 'vd';
        % disp('Software version: VD (!?)');

        % number of different scans in file stored in 2nd in
        NScans = secondInt;
        measID = fread(fid,1,'uint32');
        fileID = fread(fid,1,'uint32');
        measOffset = cell(1, NScans);
        measLength = cell(1, NScans);
        for k=1:NScans
            measOffset{k} = fread(fid,1,'uint64');
            measLength{k} = fread(fid,1,'uint64'); 
            fseek(fid, 152 - 16, 'cof');
        end
    else
        % in VB versions, the first 4 bytes indicate the beginning of the
        % raw data part of the file
        version  = 'vb';
        % disp('Software version: VB (!?)');
        measOffset{1} = 0;
        measLength{1} = fileSize;
        NScans     = 1; % VB does not support multiple scans in one file
    end
    prot = cell(1,NScans);
    for s=1:NScans
        cPos = measOffset{s};
        fseek(fid,cPos,'bof');
        hdr_len = fread(fid, 1,'uint32');
        prot{s} = read_twix_hdr(fid);
    end
    fclose(fid);
    if NScans == 1
        prot = prot{1};
    end

    
end
    