function d=b0_crct_bmir_epi(d,df,para)
% input d is in k-space
    if isvector(df)
        nx=1;
        nshot=length(df);
    else
        nx=size(df,1);
        nshot=numel(df)/nx;
    end
    if para.isgre
        % gre
        te=para.te_contr(:)*1e-3;
        dp=2*pi*te.'.*...
           reshape(df,[nx,1,1,1,nshot]);
    else
% $$$         error('*** Non-GRE data not supported yet! ***');
        % epi, echo shifting
        % 1. base line echo time for each line
        te=para.te_ro(:)-para.te+para.te_contr(:).';
        te=te(:)*1e-3;
        % 2. consider the order how echo time shift
        n_interl=para.n_interleaves;
        necho=length(para.te_contr);
        nshot=length(df);
        if isfield(para,'pe_order') && strcmp(para.pe_order,'rev_linear')
            % for philips, the pe order can be reversed
            te=te+...
               [n_interl-1:-1:0]*...
               para.echo_spacing*1e-6/n_interl*para.int_te_shift;
        else
            te=te+...
               [0:n_interl-1]*...
               para.echo_spacing*1e-6/n_interl*para.int_te_shift;
        end
        n_par_rep=para.n_partitions/para.sense_rate_s*para.n_reps;
        if isfield(para,'loop_order') && strcmp(para.loop_order,'zy_order')
            te=repmat(te,[n_par_rep,1]);
        else
            te=repmat(te,[1,n_par_rep]);
        end
        te=reshape(te,[1,para.nk_shot*necho,1,1,nshot]);
        dp=2*pi*te.*reshape(df,[1,1,1,1,nshot]);
    end
    d=fftmr(d,-1,1);
    d=d./exp(1i*dp);
end
