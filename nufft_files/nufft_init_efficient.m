function st = nufft_init_efficient(om, Nd, Jd, Kd, precision, varargin)
% Modified by Jiaen 01/2020: return sparse matrix index and value
%function st = nufft_init(om, Nd, Jd, Kd, [n_shift,] ...)
%
% Initialize structure for d-dimension NUFFT using KB interpolator,
% particularly the interpolation matrix in sparse format.
% caution: this routine can require a lot of memory!
% in
%	om [M,d]	"digital" frequencies in radians
%	Nd [d]		image dimensions (N1,N2,...,Nd)
%	Jd [d]		# of neighbors used (in each direction)
%	Kd [d]		FFT sizes (should be >= N1,N2,...)
% optional arguments
%	n_shift [d]	n = 0-n_shift to N-1-n_shift (must be first)
%	'minmax:kb'	minmax interpolator with excellent KB scaling!
%				(minmax:kb is recommended, and used by default)
%	'minmax:tuned'	minmax interpolator, somewhat numerically tuned
%	'minmax:user'	minmax interpolator with user ({alpha}, {beta})
%	'uniform'	uniform scaling factors (not recommended)
%	'kaiser'	kaiser-bessel (KB) interpolator (minmax best alpha, m)
%			or 'kaiser', alpha, m to specify parameter (vectors)
%	'linear'	linear interpolator (a terrible straw man)
%	kernel		user-provided inline interpolation kernel(k,J)
%			(or a cell array of kernels, one for each dimension)
%	'table'		use table-based interpolation rather than sparse matrix.
%			this can save a lot of memory for large problems.
%			example ..., 'table', 2^11, 'minmax:kb'
%			where 2^11 is the table over-sampling factor.
% out
%	st.p		[M, *Nd]	sparse interpolation matrix
%					(or empty if table-based) 
%	st.sn		[(Nd)]		scaling factors
%	st.Nd,Jd,Kd,om	copies of inputs
%
% *Nd is shorthand for prod(Nd).
% (Nd) is shorthand for (N1,N2,...,Nd)
%
% Like fft(), the NUFFT expects the signals to be x(0,0), ...
% Use n_shift = [N1/2, N2/2, ...] for x(-N1/2,-N2/2,...), ...
%
% Copyright 2002-5-30	Jeff Fessler	The University of Michigan

if nargin < 4, help(mfilename), error args, end

% dimensionality of input space (usually 2 or 3)
dd = length(Nd);
if dd ~= length(Jd) | dd ~= length(Kd)
	error 'inconsistent dim'
end
if dd ~= size(om,2), error(sprintf('omega needs %d columns', dd)), end
precision_double=strcmp('double',precision);

%
% process optional arguments
%

% n_shift argument? (must be first)
if length(varargin) > 0 & isnumeric(varargin{1})
	n_shift = varargin{1};
	if dd ~= length(n_shift)
		error(sprintf('n_shift needs %d columns', dd))
	end
	varargin = {varargin{2:end}};
else
	n_shift = zeros(size(Nd));
end
st.n_shift = n_shift;

% default/recommended interpolator is minmax with KB scaling factors
if length(varargin) == 0
	varargin = {'minmax:kb'};
end

st.alpha = {};
st.beta = {};
is_kaiser_scale = logical(0);

% table based?
if ischar(varargin{1}) & streq(varargin{1}, 'table')
	st = nufft_table_init(om, Nd, Jd, Kd, n_shift, varargin{2:end});
	return
end

ktype = varargin{1};

% cell array of kernel functions: {kernel1, kernel2, ..., kernelD}
if isa(ktype, 'cell')
	if isa(ktype{1}, 'inline') | isa(ktype{1}, 'function_handle')
		ktype = 'inline';
		if length(varargin) > 1, error 'excess arguments?', end
		if length(varargin{1}) ~= dd, error 'wrong # of kernels', end
		st.kernel = varargin{1};
	else
		error 'cell array should be inline kernels!?'
	end

% or a single inline kernel for all dimension
elseif isa(ktype, 'inline') | isa(ktype, 'function_handle')
	ktype = 'inline';
	if length(varargin) > 1, error 'excess arguments?', end
	for id = 1:dd
		st.kernel{id} = varargin{1};	% all same
	end

% or a string that describes the type of interpolator
elseif ~ischar(ktype)
	error 'non-string kernel type?'

end

st.ktype = ktype;

%
% set up whatever is needed for each interpolator
%
if streq(ktype, 'inline')
	% already did it above

% linear interpolator straw man
elseif streq(ktype, 'linear')
	ktype = 'inline';
	kernel = inline('(1 - abs(k/(J/2))) .* (abs(k) < J/2)', 'k', 'J');
	for id = 1:dd
		st.kernel{id} = kernel;
	end

% KB interpolator
elseif streq(ktype, 'kaiser')
	is_kaiser_scale = logical(1);

	% with minmax-optimized parameters
	if length(varargin) == 1
		for id = 1:dd
			[st.kernel{id} st.kb_alf(id) st.kb_m(id)] = ...
				kaiser_bessel('inline', Jd(id));
		end

	% with user-defined parameters
	elseif length(varargin) == 3
		alpha_list = varargin{2};
		m_list = varargin{3};
		if (length(alpha_list) ~= dd) | (length(m_list) ~= dd)
			error 'need #dim alpha and m'
		end
		for id = 1:dd
			[st.kernel{id} st.kb_alf(id) st.kb_m(id)] = ...
				kaiser_bessel('inline', Jd(id), ...
					alpha_list(id), m_list(id));
		end
	else
		error 'kaiser should have no arguments, or both alpha and m'
	end

% minmax interpolator with KB scaling factors (recommended default)
elseif streq(ktype, 'minmax:kb')
	for id = 1:dd
		[st.alpha{id}, st.beta{id}] = ...
			nufft_alpha_kb_fit(Nd(id), Jd(id), Kd(id));
	end

% minmax interpolator with numerically "tuned" scaling factors
elseif streq(ktype, 'minmax:tuned')
	for id = 1:dd
		[st.alpha{id}, st.beta{id}, ok] = ...
			nufft_best_alpha(Jd(id), 0, Kd(id)/Nd(id));
		if ~ok, error 'unknown J,K/N', end
	end

% minmax interpolator with user-provided scaling factors
elseif streq(ktype, 'minmax:user')
	if length(varargin) ~= 3, error 'user must provide alpha/beta', end
	st.alpha = varargin{2};
	st.beta = varargin{3};
	if length(st.alpha) ~= dd | length(st.beta) ~= dd
		error 'alpha/beta size mismatch'
	end

elseif streq(ktype, 'uniform')
	for id = 1:dd
		st.alpha{id} = 1;
		st.beta{id} = 0;
	end

else
	error 'unknown kernel type'
end

st.tol	= 0;

st.Jd	= Jd;
nnz_row=prod(Jd);
st.nnz_row=nnz_row;
st.Nd	= Nd;
st.Kd	= Kd;

M = size(om,1);
st.M	= M;
st.om	= om;

%
% scaling factors: "outer product" of 1D vectors
%
st.sn = 1;
for id=1:dd
	if is_kaiser_scale
		nc = [0:Nd(id)-1]'-(Nd(id)-1)/2;
		tmp = 1 ./ kaiser_bessel_ft(...
			nc/Kd(id), Jd(id), st.kb_alf(id), st.kb_m(id), 1);
	elseif streq(ktype, 'inline')
		tmp = 1 ./ nufft_interp_zn(0, Nd(id), Jd(id), Kd(id), st.kernel{id});
	else
		tmp = nufft_scale(Nd(id), Kd(id), st.alpha{id}, st.beta{id});
	end
	st.sn = st.sn(:) * tmp';
end
if length(Nd) > 1
	st.sn = reshape(st.sn, Nd);	% [(Nd)]
else
	st.sn = st.sn(:);	% [(Nd)]
end

%
% [J?,M] interpolation coefficient vectors.  will need kron of these later
%
for id=1:dd
	N = Nd(id);
	J = Jd(id);
	K = Kd(id);
	if isvar('st.kernel')
		[c, arg] = ...
		nufft_coef(om(:,id), J, K, st.kernel{id});	% [J?,M]
	else
		alpha = st.alpha{id};
		beta = st.beta{id};
		T = nufft_T(N, J, K, st.tol, alpha, beta);	% [J?,J?]
		[r, arg] = ...
		nufft_r(om(:,id), N, J, K, alpha, beta);	% [J?,M]
		if precision_double
            c = double(T * r);
        else
            c = single(T * r);
        end
        
        clear T r;
	end

	gam = 2*pi/K;
	phase_scale = 1i * gam * (N-1)/2;
    if precision_double
        phase = double(exp(phase_scale * arg));	% [J?,M] linear phase
    else
        phase = single(exp(phase_scale * arg));	% [J?,M] linear phase
    end
	
    %% test purpuse %%
% $$$     c(:)=0;
% $$$     c(ceil(J/2),:)=1;
% $$$     phase(:)=0;
% $$$     phase(ceil(J/2),:)=1;
    %% ------------ %%
	ud{id} = phase .* c;		% [J?,M]
	%
	% indices into oversampled FFT components
	%
	koff = nufft_offset(om(:,id), J, K);	% [M,1] to leftmost near nbr
	kd{id} = mod(outer_sum([1:J]', koff'), K) + 1;	% [J?,M] {1,...,K?}
	if id > 1	% trick: pre-convert these indices into offsets!
		kd{id} = (kd{id}-1) * prod(Kd(1:(id-1)));
	end

end
clear c arg gam phase phase_scale koff N J K;

%
% build sparse matrix that is [M,*Kd]
% with *Jd nonzero entries per frequency point
%

if dd >= 3
    if precision_double
        mem_float=8;
    else
        mem_float=4;
    end
     % 8 is for column and row index
    mem=nnz_row*M/1024^3*(mem_float+16);
    disp(num2str(mem,'Needs at least %g Gbyte RAM for nuFFT matrix'));
end

st.p.kk = kd{1};	% [J1,M]
st.p.uu = ud{1};	% [J1,M]
for id = 2:dd
	Jprod = prod(Jd(1:id));
	st.p.kk = block_outer_sum(st.p.kk, kd{id});	% outer sum of indices
	st.p.kk = reshape(st.p.kk, Jprod, M);
	st.p.uu = block_outer_prod(st.p.uu, ud{id});	% outer product of coefficients
	st.p.uu = reshape(st.p.uu, Jprod, M);
end	% now kk and uu are [*Jd, M]

%
% apply phase shift
% pre-do Hermitian transpose of interpolation coefficients
% [1,M]
if precision_double 
    phase = exp(1i * double((om) * n_shift(:))).';
else
    phase = exp(1i * single((om) * n_shift(:))).';
end
st.p.uu = conj(st.p.uu) .* phase(ones(1,nnz_row),:);		% [*Jd,M]

% ----------------------

st.p.mm = [1:M]; 
st.p.mm = st.p.mm(ones(nnz_row,1),:);		% [*Jd,M]
st.N=prod(Kd);
st.p.mm = reshape(st.p.mm,[nnz_row*M,1]);
st.p.kk = reshape(st.p.kk,[nnz_row*M,1]);
st.p.uu = reshape(st.p.uu,[nnz_row*M,1]);

% matlab does not allow single float
% $$$ st.p = sparse(mm(:), kk(:), uu(:), M, prod(Kd));	% sparse matrix

% in
%	x1	[J1,M]
%	x2	[J2,M]
% out
%	y	[J1,J2,M]	y(i1,i2,m) = x1(i1,m) + x2(i2,m)
%
function y = block_outer_sum(x1, x2)
[J1 M] = size(x1);
[J2 M] = size(x2);
xx1 = reshape(x1, [J1 1 M]);	% [J1,1,M] from [J1,M]
% $$$ xx1 = xx1(:,ones(J2,1),:);	% [J1,J2,M], emulating ndgrid
xx2 = reshape(x2, [1 J2 M]);	% [1,J2,M] from [J2,M]
% $$$ xx2 = xx2(ones(J1,1),:,:);	% [J1,J2,M], emulating ndgrid
y = xx1 + xx2;			% [J1,J2,M]

function y = block_outer_prod(x1, x2)
[J1 M] = size(x1);
[J2 M] = size(x2);
xx1 = reshape(x1, [J1 1 M]);	% [J1,1,M] from [J1,M]
% $$$ xx1 = xx1(:,ones(J2,1),:);	% [J1,J2,M], emulating ndgrid
xx2 = reshape(x2, [1 J2 M]);	% [1,J2,M] from [J2,M]
% $$$ xx2 = xx2(ones(J1,1),:,:);	% [J1,J2,M], emulating ndgrid
y = xx1 .* xx2;			% [J1,J2,M]
