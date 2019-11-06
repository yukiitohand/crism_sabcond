function [logYif_cor,logt_est,logAB,logBg,logIce,logYif_isnan,X,badspc]...
    = sabcondc_v3l1_gpu_batch(logYif,WA,Alib,logT,BP,varargin)
% [logYif_cor,logt_est,logAB,logBg,logIce,logYif_isnan,X,badspc]...
%     = sabcondc_v3l1_gpu_batch(logYif,WA,Alib,logT,BP,varargin)
%
% Perform sabcondc_v3l1 with Algorithm 2 (bad entries are exactly ignored).
% BATCH mode is also supported (GPU is needed). BATCH mode is activated by
% default if the third dimension of input matrices is greater than 1.
%
% INPUT PARAMETERS
%   logYif: array, [B x L x S]
%       each of the page is the observation
%   WA: array, [B x 1 x S]
%       wavelength frame
%   Alib: array, [B x Nlib x S]
%       each of the page is library matrix for column s
%   logT: array, [B x Ntc x S]
%       each of the page is collection of transmission spectra for column s
%   BP: boolean, [B x 1 x S]
%       bad pixel information
%
% OUTPUT PARAMETERS
%   logt_est: array, [B x 1 x S]
%       estimated transmission spectrum
%   logYifc_cor: array, [B x L x S]
%       corrected i/f spectra
%   logAB: array, [B x L x S]
%       estimated absorption spectra
%   logBg: array, [B x L x S]
%       estimated background spectra
%   logIce: array, [B x L x S]
%       estimated ice contributions
%   logYif_isnan: boolean array, [B x L x S]
%       bad entries are flagged.
%   X: [(1+Nlib+Nice) x L x S]
%       estimated abundance matrix
%   badspc: boolean array, [1 x L x S]
%       flag if the spectra has too many bad entries (>THREHOLD_BADSPC)
%
% OPTIONAL PARAMETERS
%  ## GENERAL PARAMETERS #-------------------------------------------------
%   'AICELIB': array, [L x Nc x S]
%       ice absporption library
%       (default) []
%   'NITER': integer
%       the number of outer iterations
%       (default) 5
%   'THRESHOLD_BADSPC': scalar
%       threshold value for which the spectra is considered to be
%       completely corrupted
%       (default) 0.8
%  ## HUWACB PARAMETERS #--------------------------------------------------
%   'LAMBDA_A': scalar,
%        trade-off parameter, controling sparsity of the coefficients of Alib
%        (default) 0.01
%   'LAMBDA_A_ICE': scalar,
%        trade-off parameter, controling sparsity of the coefficients of
%        Aicelib
%        (default) 0.0
%   'MAXITER_HUWACB': integer
%        maximum number of iteration for HUWACB
%        (default) 100
%   'TOL_HUWACB': scalar
%        tolerance value for HUWACB
%        (default) 1e-4
%   'VERBOSE_HUWACB': {'yes','no'}
%        whether or not to print the result of HUWACB or not 
%        (default) 'no'
%   'DEBUG_HUWACB': boolean
%        whether or not to debug HUWACB or not 
%        (default) false
%  ## LAD PARAMETERS #-----------------------------------------------------
%   'MAXITER_LAD': integer
%        maximum number of iteration for LAD
%        (default) 1000
%   'TOL_LAD': scalar
%        tolerance value for LAD
%        (default) 1e-4
%   'VERBOSE_LAD': {'yes','no'}
%        whether or not to print the result of LAD or not 
%        (default) 'no'
%   'DEBUG_LAD': boolean
%        whether or not to debug LAD or not 
%        (default) false
%  ## PROCESSING PARAMETERS #----------------------------------------------
%   'GPU': boolean
%        whether or not to use GPU or GPU
%        (default) true
%   'BATCH': boolean
%        whether or not to use BATCH mode processing. Only valid if 'GPU'
%        is set true.
%        (default) true if the input matrices has multiple pages, false
%                  otherwise
%   'PRECISION': {'single','doulbe'}
%        data type used for calculation
%        (default) 'double'

%% GET SIZE
[B,L,S]    = size(logYif);
[~,Nlib,~] = size(Alib);
[~,Ntc,~]  = size(logT);

%% VARARGIN
% ## GENERAL PARAMETERS #--------------------------------------------------
Aicelib   = [];
nIter     = int32(5);
th_badspc = 0.8;
% ## HUWACB PARAMETERS #---------------------------------------------------
lambda_a       = 0.01;
lambda_a_ice   = 0;
maxiter_huwacb = int32(100);
tol_huwacb     = 1e-4;
verbose_huwacb = 'no';
debug_huwacb   = false;
% ## LAD PARAMETERS #------------------------------------------------------
maxiter_lad = int32(1000);
tol_lad     = 1e-4;
verbose_lad = 'no';
debug_lad   = false;
% ## PROCESSING PARAMETERS #-----------------------------------------------
gpu       = true;
batch     = S>1;
precision = 'double';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            % ## GENERAL PARAMETERS #--------------------------------------
            case 'AICELIB'
                Aicelib = varargin{i+1};
            case 'NITER'
                nIter = varargin{i+1};
            case 'THRESHOLD_BADSPC'
                th_badspc = varargin{i+1};
            % ## HUWACB PARAMETERS #---------------------------------------
            case 'LAMBDA_A'
                lambda_a = varargin{i+1};
            case 'LAMBDA_A_ICE'
                lambda_a_ice = varargin{i+1};
            case 'MAXITER_HUWACB'
                maxiter_huwacb = round(varargin{i+1});
                if (maxiter_huwacb <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'TOL_HUWACB'
                tol_huwacb = varargin{i+1};
            case 'VERBOSE_HUWACB'
                verbose_huwacb = varargin{i+1};
            case 'DEBUG_HUWACB'
                debug_huwacb = varargin{i+1};
            % ## LAD PARAMETERS #------------------------------------------
            case 'MAXITER_LAD'
                maxiter_lad = round(varargin{i+1});
                if (maxiter_huwacb <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'TOL_LAD'
                tol_lad = varargin{i+1};
            case 'VERBOSE_LAD'
                verbose_lad = varargin{i+1};
            case 'DEBUG_LAD'
                debug_lad = varargin{i+1};
            % ## PROCESSING PARAMETERS #-----------------------------------
            case 'GPU'
                gpu = varargin{i+1};
            case 'BATCH'
                batch = varargin{i+1};
            case 'PRECISION'
                precision = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end


%%
%-------------------------------------------------------------------------%
% Preprocessing
%-------------------------------------------------------------------------%
[~,Nice,~] = size(Aicelib);

if batch
    gpu_varargin = {'gpuArray'};
else
    gpu_varargin = {};
end
if batch && ~gpu, error('BATCH Processing only works when GPU is set true'); end
if batch && S==1, fprintf('BATCH Processing is efficient if you have 3D array.'); end

if Nice==0
    Aicelib = zeros(B,0,S,precision,gpu_varargin{:});
end

% preprocessing is necessary because some of the values in logYif might be
% nans due to harse noise.
logYif_isnan = isnan(logYif);
% replace those values with linear interpolation
% logYif = interp_nan_column(logYif,logYif_isnan,WA);
% mark bad spectra
badspc = (sum(logYif_isnan,1)/B) > 0.5;

% mark bad pixels and badspc
logYif_isnan = or(logYif_isnan,BP);
logYif_isnan = or(logYif_isnan,badspc);

% save the original bad entries for future use.
logYif_isnan_ori = logYif_isnan;
logYif(logYif_isnan) = 0;
% tic; 
for i=1:S
    logYif(:,:,i) = interp_nan_column(logYif(:,:,i),logYif_isnan(:,:,i),WA(:,i));
end
% toc;

if batch
    logYif  = gpuArray(logYif);
    WA      = gpuArray(WA);
    Alib    = gpuArray(Alib);
    Aicelib = gpuArray(Aicelib);
    logT    = gpuArray(logT);
    BP      = gpuArray(BP);
    logYif_isnan = gpuArray(logYif_isnan);
    logYif_isnan_ori = gpuArray(logYif_isnan_ori);
    badspc = gpuArray(badspc);
end

%%
%-------------------------------------------------------------------------%
% initialization using the matrix logT
%-------------------------------------------------------------------------%
A = cat(2,logT,Aicelib,Alib);
% compute concave basese beforehand
% C = concaveOperator(WA);
% Dinv = C \ eye(B);
if batch
    C = concaveOperator_v2(WA);
    Dinv = pagefun(@mldivide,C,eye(B,gpu_varargin{:}));
    s_d = vnorms(Dinv,1);
    C = Dinv./s_d;
    clear Dinv s_d;
else
    C = continuumDictionary(WA);
    s_c = vnorms(C,1);
    C = C./s_c;
    C = C*2;
    clear s_c;
end
[~,Nc,~] = size(C);
if strcmpi(precision,'single')
    C = single(C);
end

% create lambdas
% lambda_c = zeros(B,L,S,precision);
% lambda_a2 = zeros(Nlib+Ntc,L,S,precision);
% lambda_r = ones(B,L,S,precision);
lambda_c = zeros(B,L,S,precision,gpu_varargin{:});
lambda_a_2 = zeros(Nlib+Nice+Ntc,L,S,precision,gpu_varargin{:});
lambda_r = ones(B,L,S,precision,gpu_varargin{:});

lambda_a_2((1+Ntc):(Ntc+Nice),:,:) = lambda_a_ice.*ones(Nice,L,S,precision,gpu_varargin{:});
lambda_a_2((1+Ntc+Nice):end,:,:) = lambda_a.*ones(Nlib,L,S,precision,gpu_varargin{:});

% 
lambda_c(logYif_isnan) = inf;
lambda_c([1,Nc],:,:) = 0; % safeguard
lambda_r(logYif_isnan) = 0;

c2_z = zeros([Nc,1],precision);
% c2_z = zeros([Nc,1],precision,'gpuArray');
c2_z(1) = -inf; c2_z(Nc) = -inf;

% main computation
% tic;
if batch
    [ X,Z,~,D,rho,Rhov,~,~,cost_val]...
        = huwacbl1_admm_gat_a_batch(A,logYif,C,'LAMBDA_A',lambda_a_2,...
                'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,'C2_Z',c2_z,......
                'tol',1e-5,'maxiter',maxiter_huwacb,'verbose',verbose_huwacb,...
                'precision',precision); % toc;
else
    [  X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
    = huwacbl1_admm_gat_a(A,logYif,WA,'LAMBDA_A',lambda_a_2,...
    'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
    'tol',1e-5,'maxiter',maxiter_huwacb,'verbose','no',...
    'precision',precision,'gpu',gpu,'Concavebase',C,'debug',debug_huwacb);
end
% toc;

% evaluate bad pixels
if batch
    RR = logYif - pagefun(@mtimes,A,X) - pagefun(@mtimes,C,Z);
else
    RR = logYif - A*X - C*Z;
end
RR_std = nanstd(RR,[],2);
logYif_isnan = or( logYif_isnan_ori, abs(RR)>0.1 );
logYif_isnan = or(logYif_isnan,RR_std>0.015);
% finally flag spectra that have too many bad channels.
badspc = (sum(logYif_isnan,1)/B) > th_badspc;
logYif_isnan = or(logYif_isnan,badspc);


lambda_c(logYif_isnan) = inf; lambda_c(~logYif_isnan) = 0;
lambda_c([1,Nc],:,:) = 0; % safeguard
lambda_r(logYif_isnan) = 0; lambda_r(~logYif_isnan) = 1;

resNrm = nansum(abs(RR).* lambda_r,[1,2]);

% Get initial transmission spectrum
if batch
    RR  = RR + pagefun(@mtimes,logT,X(1:Ntc,:,:));
else
    RR = RR + logT*X(1:Ntc,:,:);
end
Xtc = sum(X(1:Ntc,:,:),1);
if batch
    [logt_est,r_lad,d_lad,rho_lad,Rhov_lad,~,~,cost_val,Kcond]...
           = lad_admm_gat_b_batch(permute(Xtc,[2,1,3]), permute(RR,[2,1,3]),...
                'lambda_r',permute(lambda_r,[2,1,3]),'tol',tol_lad,'maxiter',maxiter_lad,...
                'verbose',verbose_lad,'precision',precision);
else
    [logt_est,~,~,rho_lad,Rhov_lad,~,~,cost_val]...
    = lad_admm_gat_b(permute(Xtc,[2,1,3]), permute(RR,[2,1,3]),...
             'lambda_r',permute(lambda_r,[2,1,3]),...
             'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad,...
             'PRECISION',precision,'gpu',gpu,'debug',debug_lad);
end
%
logt_est = permute(logt_est,[2,1,3]);
if batch
    RR = RR - pagefun(@mtimes,logt_est,Xtc);
else
    RR = RR - logt_est*Xtc;
end
resNewNrm = nansum(abs(lambda_r .* RR),[1,2]);

%%
%-------------------------------------------------------------------------%
% main loop
%-------------------------------------------------------------------------%
A = cat(2,logt_est,Alib);
X = cat(1,Xtc,X((1+Ntc):end,:,:));
D = cat(1,zeros(1,L,S,precision,gpu_varargin{:}),D((1+Ntc):end,:,:));
lambda_a_2 = ones(1+Nice+Nlib,L,S,precision,gpu_varargin{:});
lambda_a_2(1,:,:) = 0;
lambda_a_2(2:(Nice+1),:,:) = lambda_a_ice.*ones(Nice,L,S,precision,gpu_varargin{:});
lambda_a_2((2+Nice):end,:,:) = lambda_a.*ones(Nlib,L,S,precision,gpu_varargin{:});
% rho = ones([1,L,S],precision,'gpuArray');
Rhov = cat(1,ones(1,1,S,precision,gpu_varargin{:}),Rhov((Ntc+1):(Ntc+Nlib+Nc+B),:,:));
% lambda_tmp = lambda_a;
% always update lambda_tmp
% lambda_tmp = lambda_tmp .* resNewNrm ./ resNrm;

for n=2:nIter
    lambda_a_2(2:end,:,:) = lambda_a_2(2:end,:,:) .* resNewNrm ./ resNrm;
    % tic;
    if n==2
        if batch
            [ X,Z,~,D,rho,Rhov,~,~,cost_val,Tcond ]...
                = huwacbl1_admm_gat_a_batch(A,logYif,C,...
                    'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
                    'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
                    'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb,...
                    'PRECISION',precision);
        else
            [ X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
            = huwacbl1_admm_gat_a(A,logYif,WA,'LAMBDA_A',lambda_a_2,...
            'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
            'Z0',Z,'D0',D,'X0',X,'R0',RR,...
            'tol',1e-5,'maxiter',maxiter_huwacb,'verbose','no',...
            'precision',precision,'gpu',gpu,'Concavebase',C,...
            'debug',debug_huwacb);
        end
    else
        if batch
           [ X,Z,~,D,rho,Rhov,~,~,cost_val ] = huwacbl1_admm_gat_a_batch(A,logYif,C,...
              'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
              'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
              'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
              'PRECISION',precision,'Tcond',Tcond);
        else
            [ X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
            = huwacbl1_admm_gat_a(A,logYif,WA,'LAMBDA_A',lambda_a_2,...
            'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
            'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
            'tol',tol_huwacb,'maxiter',maxiter_huwacb,'verbose','no',...
            'precision',precision,'gpu',gpu,'Concavebase',C,...
            'debug',debug_huwacb);
        end
    end
    % toc;
    % evaluate bad pixels
    if batch
        RR = logYif - pagefun(@mtimes,A,X) - pagefun(@mtimes,C,Z);
    else
        RR = logYif - A*X - C*Z;
    end
    logYif_isnan = or( logYif_isnan_ori, abs(RR)>0.015 );
    badspc = (sum(logYif_isnan,1)/B) > th_badspc;
    logYif_isnan = or(logYif_isnan,badspc);

    lambda_c(logYif_isnan) = inf; lambda_c(~logYif_isnan) = 0;
    lambda_c([1,Nc],:,:) = 0; % safeguard
    lambda_r(logYif_isnan) = 0; lambda_r(~logYif_isnan) = 1;
    
    resNrm = nansum(abs(lambda_r .* RR),[1,2]);
    
    % update logt_est
    if batch
        RR  = RR + pagefun(@mtimes,A(:,1,:),X(1,:,:));
    else
        RR  = RR + A(:,1)*X(1,:);
    end
    if batch
        [logt_est,r_lad,d_lad,rho_lad,Rhov_lad,~,~,cost_val]...
           = lad_admm_gat_b_batch(permute(Xtc,[2,1,3]), permute(RR,[2,1,3]),...% 'rho',rho_lad,'Rhov',Rhov_lad,...
                'lambda_r',permute(lambda_r,[2,1,3]),'tol',tol_lad,'maxiter',maxiter_lad,...
                'verbose',verbose_lad,'precision',precision,'Kcond',Kcond);
    else
        [logt_est,~,~,rho_lad,Rhov_lad,~,~,cost_val]...
            = lad_admm_gat_b(permute(Xtc,[2,1,3]), permute(RR,[2,1,3]),...
             'lambda_r',permute(lambda_r,[2,1,3]),...
             'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad,...
             'PRECISION',precision,'gpu',gpu,'debug',debug_lad);
        
    end
    logt_est = permute(logt_est,[2,1,3]);
    if batch
        RR = RR - pagefun(@mtimes,logt_est,X(1,:,:));
    else
        RR = RR - logt_est*X(1,:);
    end
    resNewNrm = nansum(abs(lambda_r .* RR),[1,2]);
    
    A(:,1,:) = logt_est;
    
end

%%
%-------------------------------------------------------------------------%
% last iteration
%-------------------------------------------------------------------------%
if batch
    [ X,Z,D,rho,Rhov,~,~,cost_val ] = huwacbl1_admm_gat_a_batch(A,...
        logYif,C,...
        'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
        'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,...%'rho',rho,'Rhov',Rhov,...
        'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
        'PRECISION',precision);
else
    [ X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
            = huwacbl1_admm_gat_a(A,logYif,WA,'LAMBDA_A',lambda_a_2,...
            'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
            'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
            'tol',tol_huwacb,'maxiter',maxiter_huwacb,'verbose','no',...
            'precision',precision,'gpu',gpu,'Concavebase',C,...
            'debug',debug_huwacb);
end

% substituting all the variables
if batch
    logAB        = pagefun(@mtimes,Alib,X((2+Nice):end,:,:));
    logBg        = pagefun(@mtimes,C,Z);
    logIce       = pagefun(@mtimes,Aicelib,X(2:(Nice+1),:,:));
    logYif_cor   = logYif - pagefun(@mtimes,logt_est,X(1,:,:)) - logIce;
else
    logAB = Alib*X((2+Nice):end,:);
    logBg = C*Z;
    logIce = Aicelib*X(2:(Nice+1),:);
    logYif_cor = logYif - logt_est*X(1,:) - logIce;
end
logt_est     = A(:,1,:);
logYif_cor(logYif_isnan) = nan;

if batch
    [logYif_cor,logt_est,logAB,logBg,logIce,logYif_isnan,X,badspc]...
        = gather(logYif_cor,logt_est,logAB,logBg,logIce,logYif_isnan,X,badspc);
else
end

end

