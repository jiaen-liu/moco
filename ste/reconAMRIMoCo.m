function reconAMRIMoCo(folder,mid,fPar,steps,varargin)
% example reconAMRIMoCo ~/data/20200902_1 125 {'~/matlab/ste/steParMoCoB0Co.conf'}
    if ischar(mid)
        mid=eval(mid);
    end
    if ischar(fPar)
        fPar=eval(fPar);
    end
    if nargin>3
        if ischar(steps)
            steps=eval(steps);
        end
    else
        steps=[];
    end
    % parse input variables
    p=inputParser;
    parpool_path='~/matlab/';
    par_overwrite=[];
    addParameter(p,'parpool_path',parpool_path,@ischar);
    addParameter(p,'par_overwrite',par_overwrite,@isstruct);
    p.parse(varargin{:});
    parpool_path=p.Results.parpool_path;
    par_overwrite=p.Results.par_overwrite;
    % 
    nrecon=length(fPar);
    cd(folder);
    nmid=length(mid);
    % set up parallel pool configuration
    parpool_file=fullfile(parpool_path,'local.settings');
    if isdeployed() && exist(parpool_file)
        setmcruserdata('ParallelProfile',parpool_file);
        distcomp.feature('LocalUseMpiexec', false);
    else
        warning('*** Parpool setting file doesn''t exist! ***');
    end
    for i=1:nmid
        disp(['*** MID',num2str(mid(i)),' ***']);
        clear pari pario;
        for j=1:nrecon
            [~,fn]=fileparts(pwd);
            if exist('pari')
                pario=pari;
            end
            % retrieve parameters from the
            % configuration file
            pari=readFilePar(fPar{j});
            if ~isempty(par_overwrite)
                pari=pass_var_struct(pari,par_overwrite);
            end
            if ~isempty(steps)
                strStep=steps;
            else
                if j==1
                    % Run through all if it's the first reconstruction
                    strStep={'full'};
                else
                    % if the cluster number changes
                    % run the cluster step again
                    if pari.nc_max==pario.nc_max
                        strStep={'grerec'};
                    else
                        strStep={'cluster','grerec'};
                    end
                end
            end
            % reconstruction
            [im_recon,par,para_mr]=...
                recon_epi_ste_beta(mid(i),...
                                   pari,strStep);
            if isempty(im_recon)
                return;
            end
            % save result
            save_mat(rp(['im',...
                         num2str(mid(i)),...
                         '.result_',...
                         par.recon_type,'.mat']),'im_recon','par',...
                     'overwrite',1,'read4all',1);
            path_output=rp('');
            prefix=['im',num2str(mid(i)),'_',...
                    par.recon_type,];
            save_amri_epi_nii(path_output,im_recon,para_mr,prefix);
        end
    end
end