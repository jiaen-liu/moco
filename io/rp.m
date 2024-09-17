function [ fullname ] = rp(filename,varargin)

% 2022-10-25 PvG & JAdZ
%     Global variable amriMoCo_dirout will override the default output directory
%     if it is defined.

    nvarar=length(varargin);
    % define the parent folder
    if nvarar > 0
        parent=varargin{1};
    end
    archive=0;
    if nvarar > 1
        archive=varargin{2};
    end
    if exist('parent') 
        name=parent;
    else
        [filepath,name]=fileparts(pwd);
    end
	global amriMoCo_dirout
	if (class(amriMoCo_dirout) == "char")
		dout=amriMoCo_dirout;
	else
		dout='result';
	end
    raid_common=fullfile('/','raid','common');
    data_folder=fullfile(raid_common,name);
    % Archiving
    if archive
        if length(name)>=4
            if contains(name(1:2),'20')
                data_year=str2num(name(1:4));
                % try a few locations
                % /raid/common/byyear/year/
                if ~exist(data_folder)
                    data_folder=fullfile(raid_common,...
                                         'byyear',...
                                         name(1:4),...
                                         name);
                end
                % /raid/common/fmrif7t
                if ~exist(data_folder)
                    data_folder=fullfile(raid_common,'fmrif7t',...
                                         name);
                end
                % /raid/common/fmrif3td
                if ~exist(data_folder)
                    data_folder=fullfile(raid_common,'fmrif3td',...
                                         name);
                end
                if ~exist(data_folder)
                    error('*** The data doesn''t exist! ***');
                end
                rdir=fullfile(data_folder,'result_JL');
            else
                rdir=fullfile('~',dout,name);
            end
        else
            rdir=fullfile('~',dout,name);
        end
    else
        if contains(getenv('HOSTNAME'),'biowulf','IgnoreCase',true) || ...
            contains(getenv('HOSTNAME'),'cn','IgnoreCase',true)
            % on biowulf, /lscratch/jobid for storing temporary files
            rdir=fullfile(getenv('TMPDIR'),dout,name);
        else
            % on others, save data in the data folder with a different name
            rdir=fullfile('.',dout);
        end
    end
    if exist(rdir) ~= 7
        status=mkdir(rdir);
        if ~status
            error(['*** Unable to create directory ''',...
                   rdir,...
                   '''! ***']);
        end
        file_permission(rdir,'+rw','ugo');
    end
    fullname=[rdir filesep filename];
end

