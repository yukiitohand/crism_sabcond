function [logYif_cor,logt_est,logAB,logBg,logIce,logYif_isnan,Xt,Xlib,Xice,badspc]...
    = sabcondc_v3l1_gpu_batch(logYif,WA,Alib,logT,BP,varargin)
% [logYif_cor,logt_est,logAB,logBg,logIce,logYif_isnan,Xt,Xlib,Xice,badspc]...
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
%   BP: boolean or 1nan form, [B x 1 x S]
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
%   Xt: [1 x L x S]
%       estimated path length matrix
%   Xlib: [Nlib x L x S]
%       estimated abundance matrix associated with lib
%   Xice: [Nice x L x S]
%       estimated coefficients associated with Aicelib
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
%   'LAMBDA_UPDATE_RULE': string {'L1SUM','MED','NONE'}
%       define how to update trade-off parameters.
%       (default) 'L1SUM'
%   'FFC_MODE': boolean,
%       if true, path length is fixed to one.
%       (default) false
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
%  ## WEIGHT PARAMETERS #--------------------------------------------------
%   'WEIGHT_MODE': integer {0,1}
%       mode id of weight option. 0 for uniform weight, 1 for weight based
%       on noise estimation.
%       (default) 0
%   ********PARAMETERS BELOW ARE ONLY REQUIRED FOR 'WEIGHT_MODE'=0*********
%   'STDL1_IFDF': array, [ B x 1 x S]
%       median absolute deviation for median calculated from DF.
%       (default) []
%   'SFIMG': array, [ B x 1 x S]
%       sunflux image from CDR
%       (default) []
%   'WA_um_pitch': array, [ B x 1 x S]
%       wavelength pitch for each channel in micron unit.
%       (default) []
%   'LBL': struct
%       CRISM pds lbl struct
%       (default) []

%% GET SIZE
[B,L,S]    = size(logYif);
[~,Nlib,~] = size(Alib);
[~,Ntc,~]  = size(logT);

%% VARARGIN
% ## GENERAL PARAMETERS #--------------------------------------------------
Aicelib   = [];
nIter     = int32(5);
th_badspc = 0.8;
lambda_update_rule = 'L1SUM';
ffc_mode  = false;
% ## WEIGHT PARAMETERS #---------------------------------------------------
weight_mode = 0;
stdl1_ifdf  = [];
SFimg       = [];
WA_um_pitch = [];
lbl         = [];
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
            case 'LAMBDA_UPDATE_RULE'
                lambda_update_rule = varargin{i+1};
            case 'FFC_MODE'
                ffc_mode = varargin{i+1};
                
            % ## WEIGHT PARAMETERS #---------------------------------------
            case 'WEIGHT_MODE'
                weight_mode = varargin{i+1};
            case 'STDL1_IFDF'
                stdl1_ifdf = varargin{i+1};
            case 'SFIMG'
                SFimg = varargin{i+1};
            case 'WA_UM_PITCH'
                WA_um_pitch = varargin{i+1};
            case 'LBL'
                lbl = varargin{i+1};
            
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


switch weight_mode
    case 0
        y_normalize = true;
    case 1
        y_normalize = false;
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
BP = (BP==1);
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
    stdl1_ifdf   = gpuArray(stdl1_ifdf);
    SFimg       = gpuArray(SFimg);
    WA_um_pitch = gpuArray(WA_um_pitch);
end
Yif = exp(logYif);

%%
%-------------------------------------------------------------------------%
% initialization using the matrix logT
%-------------------------------------------------------------------------%
if ffc_mode
    A = cat(2,Aicelib,Alib);
    N_A = Nice+Nlib;
    idxAlogT = false(1,N_A);
    idxAice = false(1,N_A); idxAice(1:Nice) = true;
    idxAlib = false(1,N_A1); idxAlib((Nice+1):N_A) = true;
else
    A = cat(2,logT,Aicelib,Alib);
    N_A = Ntc+Nice+Nlib;
    idxAlogT = false(1,N_A); idxAlogT(1:Ntc) = true;
    idxAice = false(1,N_A);  idxAice((Ntc+1):(Ntc+Nice)) = true;
    idxAlib = false(1,N_A);  idxAlib((Ntc+Nice+1):N_A) = true;
end
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
if ffc_mode
    lambda_a_2 = zeros(Nice+Nlib,L,S,precision,gpu_varargin{:});
else
    lambda_a_2 = zeros(Ntc+Nice+Nlib,L,S,precision,gpu_varargin{:});
end

switch weight_mode
    case 0
        lambda_r = ones(B,L,S,precision,gpu_varargin{:});
        lambda_r(logYif_isnan_ori) = 0;
    case 1
        % compute weight
        mYif = nanmean(Yif,2);
        if batch
            mYif1 = pagefun(@mtimes,pagefun(@transpose,mYif),Yif)./ sum(mYif.^2,1);
            Ymdl = pagefun(@mtimes,mYif,mYif1);
        else
            mYif1 = mYif'*Yif/norm(mYif,2)^2;
            Ymdl = mYif*mYif1;
        end
        RDimg = if2rd(Ymdl,SFimg,lbl);
        [photon_noise_mad_stdif] = estimate_photon_noise_CRISM_base(...
                                        RDimg,WA,WA_um_pitch,lbl,SFimg);
        %
        lambda_r = 1./(stdl1_ifdf+photon_noise_mad_stdif).*(Ymdl)/(B*20);
        lambda_r(logYif_isnan_ori) = 0;
end

lambda_a_2(idxAice,:,:) = lambda_a_ice.*ones(Nice,L,S,precision,gpu_varargin{:});
lambda_a_2(idxAlib,:,:) = lambda_a.*ones(Nlib,L,S,precision,gpu_varargin{:});

% 
lambda_c(logYif_isnan) = inf;
lambda_c([1,Nc],:,:) = 0; % safeguard


c2_z = zeros([Nc,1],precision);
% c2_z = zeros([Nc,1],precision,'gpuArray');
c2_z(1) = -inf; c2_z(Nc) = -inf;

% main computation
% tic;
if ffc_mode
    if batch
        [ X,Z,~,D,rho,Rhov,~,~,cost_val]...
            = huwacbl1_admm_gat_a_batch(A,logYif-logT,C,'LAMBDA_A',lambda_a_2,...
                    'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,'C2_Z',c2_z,......
                    'tol',1e-5,'maxiter',maxiter_huwacb,'verbose',verbose_huwacb,...
                    'precision',precision,'YNORMALIZE',y_normalize);
    else
        [  X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
            = huwacbl1_admm_gat_a(A,logYif-logT,WA,'LAMBDA_A',lambda_a_2,...
            'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,'YNORMALIZE',y_normalize,.....
            'tol',1e-5,'maxiter',maxiter_huwacb,'verbose','no',...
            'precision',precision,'gpu',gpu,'Concavebase',C,'debug',debug_huwacb);
    end
else
    if batch
        [ X,Z,~,D,rho,Rhov,~,~,cost_val]...
            = huwacbl1_admm_gat_a_batch(A,logYif,C,'LAMBDA_A',lambda_a_2,...
                    'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,'C2_Z',c2_z,......
                    'tol',1e-5,'maxiter',maxiter_huwacb,'verbose',verbose_huwacb,...
                    'precision',precision,'YNORMALIZE',y_normalize); % toc;
    else
        [  X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
        = huwacbl1_admm_gat_a(A,logYif,WA,'LAMBDA_A',lambda_a_2,...
        'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,'YNORMALIZE',y_normalize,.....
        'tol',1e-5,'maxiter',maxiter_huwacb,'verbose','no',...
        'precision',precision,'gpu',gpu,'Concavebase',C,'debug',debug_huwacb);
    end
end
% toc;

% evaluate bad pixels
if batch
    % Ymdl = pagefun(@mtimes,A,X) + pagefun(@mtimes,C,Z);
    RR   = logYif - pagefun(@mtimes,A,X) - pagefun(@mtimes,C,Z);
else
    % Ymdl = A*X + C*Z;
    RR = logYif - A*X - C*Z;
end

if ffc_mode
    Ymdl = Ymdl + logT;
end

% RR = logYif - Ymdl;

% ## denoising ##----------------------------------------------------------
switch weight_mode
    case 0
        RR_std = nanstd(RR,[],2);
        % noticed this de-noising is slightly different from the one in
        % sabcondc_v3l1_pub.m
        % I think the below one is the one I wanted to try...
        % logYif_isnan = or(logYif_isnan_ori,and(abs(RR)>0.1,RR_std>0.015));
        %
        logYif_isnan = or( logYif_isnan_ori, abs(RR)>0.1 );
        logYif_isnan = or(logYif_isnan,RR_std>0.015);
        
        % finally flag spectra that have too many bad channels.
        badspc = (sum(logYif_isnan,1)/B) > th_badspc;
        logYif_isnan = or(logYif_isnan,badspc);
        lambda_r(logYif_isnan) = 0; lambda_r(~logYif_isnan) = 1;
        
    case 1
        Ymdl = exp(Ymdl);
        RDimg = if2rd(Ymdl,SFimg,lbl);
        [photon_noise_mad_stdif] = estimate_photon_noise_CRISM_base(...
                                        RDimg,WA,WA_um_pitch,lbl,SFimg);
        mad_rr_band_theor = nanmedian(stdl1_ifdf+photon_noise_mad_stdif,2);
        res_exp = Yif - Ymdl;
        mad_rr_band_prac = robust_v3('med_abs_dev_from_med',res_exp,2,'NOutliers',10);
        mad_log_band = robust_v3('med_abs_dev_from_med',RR,2,'NOutliers',10);
        mad_rr_band = max(mad_rr_band_theor,mad_rr_band_prac);
        
        bp_est_bool = and(mad_rr_band>0.0015, mad_log_band>0.005);
        
        logYif_isnan = or(logYif_isnan_ori,bp_est_bool);
        
        % check too many bad entries are detected for each spectrum.
        badspc = (sum(logYif_isnan,1)/B) > th_badspc;
        
        
        logYif_isnan = or(logYif_isnan,badspc);

        
        % update lambda_r? not sure.
        
end

lambda_c(logYif_isnan) = inf; lambda_c(~logYif_isnan) = 0;
lambda_c([1,Nc],:,:) = 0; % safeguard

%--------------------------------------------------------------------------


resNrm = nansum(abs(RR).* lambda_r,[1,2]);

% Get initial transmission spectrum
if ffc_mode
    RR = RR + logT;
    Xtc = ones(1,L,S,precision,gpu_varargin);
else
    if batch
        RR  = RR + pagefun(@mtimes,logT,X(idxAlogT,:,:));
    else
        RR = RR + logT*X(idxAlogT,:,:);
    end
    Xtc = sum(X(idxAlogT,:,:),1);
end

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

switch weight_mode
    case 0

    case 1
        if batch
            Ymdl = exp(log(Ymdl) - pagefun(@mtimes,logT,X(idxlogT,:,:)) ...
                + pagefun(@mtimes,logt_est,Xtc));
        else
            Ymdl = exp(log(Ymdl) - logT*X(idxlogT,:,:) + logt_est*Xtc);
        end
        RDimg = if2rd(Ymdl,SFimg,lbl);
        [photon_noise_mad_stdif] = estimate_photon_noise_CRISM_base(...
                                        RDimg,WA,WA_um_pitch,lbl,SFimg);
        mad_rr_theor = stdl1_ifdf+photon_noise_mad_stdif;
        res_exp = Yif - Ymdl;
        mad_rr_band_prac = robust_v3('med_abs_dev_from_med',res_exp,2,...
            'NOutliers',10);
        
        lambda_r = 1./max(mad_rr_theor,mad_rr_band_prac).*(Ymdl)./(B*20);
        
        % do we need this?
        lambda_r(logYif_isnan) = 0;
end

resNewNrm = nansum(abs(lambda_r .* RR),[1,2]);

%%
%-------------------------------------------------------------------------%
% main loop
%-------------------------------------------------------------------------%
if ffc_mode
    A = cat(2,Aicelib,Alib);
    % X = X; D = D;
    N_A = Nice+Nlib;
    idxAlogT = false(1,N_A);
    idxAice = false(1,N_A);  idxAice(1:Nice) = true;
    idxAlib = false(1,N_A);  idxAlib((Nice+1):N_A) = true;
    lambda_a_2 = ones(N_A,L,S,precision,gpu_varargin{:});
    lambda_a_2(idxAice,:,:) = lambda_a_ice.*ones(Nice,L,S,precision,gpu_varargin{:});
    lambda_a_2(idxAlib,:,:) = lambda_a.*ones(Nlib,L,S,precision,gpu_varargin{:});
else
    A = cat(2,logt_est,Aicelib,Alib);
    X = cat(1,Xtc,X((1+Ntc):N_A,:,:));
    D = cat(1,zeros(1,L,S,precision,gpu_varargin{:}),D((1+Ntc):(N_A+Nc+B),:,:));
    Rhov = cat(1,ones(1,1,S,precision,gpu_varargin{:}),Rhov((Ntc+1):(N_A+Nc+B),:,:));
    
    N_A = 1+Nice+Nlib;
    idxAlogT = false(1,N_A); idxAlogT(1) = true;
    idxAice = false(1,N_A);  idxAice(2:(Nice+1)) = true;
    idxAlib = false(1,N_A);  idxAlib((Nice+2):N_A) = true;
    lambda_a_2 = ones(N_A,L,S,precision,gpu_varargin{:});
    lambda_a_2(idxAlogT,:,:) = 0;
    lambda_a_2(idxAice,:,:) = lambda_a_ice.*ones(Nice,L,S,precision,gpu_varargin{:});
    lambda_a_2(idxAlib,:,:) = lambda_a.*ones(Nlib,L,S,precision,gpu_varargin{:});
end

switch upper(lambda_update_rule)
    case 'L1SUM'
        lambda_a_2(2:end,:,:) = lambda_a_2(2:end,:,:) .* resNewNrm ./ resNrm;
    case 'MED'
        error('not implemented yet');
    case 'NONE'
    otherwise
        error('Undefined LAMBDA_UPDATE_RULE: %s',lambda_update_rule);
end
% rho = ones([1,L,S],precision,'gpuArray');

% lambda_tmp = lambda_a;
% always update lambda_tmp
% lambda_tmp = lambda_tmp .* resNewNrm ./ resNrm;

for n=2:nIter
    % tic;
    if n==2
        if ffc_mode
            if batch
                 [ X,Z,~,D,rho,Rhov,~,~,cost_val,Tcond ]...
                    = huwacbl1_admm_gat_a_batch(A,logYif-logt_est,C,...
                        'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
                        'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
                        'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb,...
                        'PRECISION',precision,'YNORMALIZE',y_normalize);
            else
                [ X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
                = huwacbl1_admm_gat_a(A,logYif-logt_est,WA,'LAMBDA_A',lambda_a_2,...
                'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
                'Z0',Z,'D0',D,'X0',X,'R0',RR,...
                'tol',1e-5,'maxiter',maxiter_huwacb,'verbose','no',...
                'precision',precision,'gpu',gpu,'Concavebase',C,...
                'debug',debug_huwacb,'YNORMALIZE',y_normalize);
            end
        else
            if batch
                [ X,Z,~,D,rho,Rhov,~,~,cost_val,Tcond ]...
                    = huwacbl1_admm_gat_a_batch(A,logYif,C,...
                        'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
                        'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
                        'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb,...
                        'PRECISION',precision,'YNORMALIZE',y_normalize);
            else
                [ X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
                = huwacbl1_admm_gat_a(A,logYif,WA,'LAMBDA_A',lambda_a_2,...
                'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
                'Z0',Z,'D0',D,'X0',X,'R0',RR,...
                'tol',1e-5,'maxiter',maxiter_huwacb,'verbose','no',...
                'precision',precision,'gpu',gpu,'Concavebase',C,...
                'debug',debug_huwacb,'YNORMALIZE',y_normalize);
            end
        end
    else
        if ffc_mode
            if batch
                [ X,Z,~,D,rho,Rhov,~,~,cost_val ] = huwacbl1_admm_gat_a_batch(A,logYif-logt_est,C,...
                  'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
                  'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
                  'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
                  'PRECISION',precision,'Tcond',Tcond,'YNORMALIZE',y_normalize);
            else
                [ X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
                = huwacbl1_admm_gat_a(A,logYif-logt_est,WA,'LAMBDA_A',lambda_a_2,...
                'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
                'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
                'tol',tol_huwacb,'maxiter',maxiter_huwacb,'verbose','no',...
                'precision',precision,'gpu',gpu,'Concavebase',C,...
                'debug',debug_huwacb,'YNORMALIZE',y_normalize);
            end
        else
            if batch
               [ X,Z,~,D,rho,Rhov,~,~,cost_val ] = huwacbl1_admm_gat_a_batch(A,logYif,C,...
                  'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
                  'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
                  'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
                  'PRECISION',precision,'Tcond',Tcond,'YNORMALIZE',y_normalize);
            else
                [ X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
                = huwacbl1_admm_gat_a(A,logYif,WA,'LAMBDA_A',lambda_a_2,...
                'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
                'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
                'tol',tol_huwacb,'maxiter',maxiter_huwacb,'verbose','no',...
                'precision',precision,'gpu',gpu,'Concavebase',C,...
                'debug',debug_huwacb,'YNORMALIZE',y_normalize);
            end
        end
    end
    % toc;
    % evaluate bad pixels
    if batch
        % Ymdl = pagefun(@mtimes,A,X) + pagefun(@mtimes,C,Z);
        RR = logYif - pagefun(@mtimes,A,X) - pagefun(@mtimes,C,Z);
    else
        % Ymdl = A*X + C*Z;
        RR = logYif - A*X - C*Z;
    end
    if ffc_mode
        Ymdl = Ymdl + logT;
    end
    %RR = logYif - Ymdl;
    
    % ## denoising ##------------------------------------------------------
    switch weight_mode
        case 0
            logYif_isnan = or( logYif_isnan_ori, abs(RR)>0.015 );
            badspc = (sum(logYif_isnan,1)/B) > th_badspc;
            logYif_isnan = or(logYif_isnan,badspc);
            lambda_r(logYif_isnan) = 0; lambda_r(~logYif_isnan) = 1;
        case 1
            Ymdl = exp(Ymdl);
            RDimg = if2rd(Ymdl,SFimg,lbl);
            [photon_noise_mad_stdif] = estimate_photon_noise_CRISM_base(...
                                        RDimg,WA,WA_um_pitch,lbl,SFimg);
            res_exp = Yif - Ymdl;
            mad_rr_band_prac = robust_v3('med_abs_dev_from_med',res_exp,2,'NOutliers',10);
            mad_rr_band = max(mad_rr_band_theor,mad_rr_band_prac);
            mad_expected = max(mad_rr_band_prac,(stdl1_ifdf+photon_noise_mad_stdif));

            % More rigid bad pixel detection
            bp_est_bool = mad_rr_band>0.001;
            
            % Now perform temporal spike removal (assuming that bias 
            % problem is gone)
            % The residual is evaluated in atm-corrected i/f domain.
            res_exp_scaled = res_exp./(exp(A(:,1)).^X(1,:));
            logYif_isnan_spk = abs(res_exp_scaled)>0.0015;
            
            % combine bad entry detections
            logYif_isnan = or(logYif_isnan_ori, or(bp_est_bool,logYif_isnan_spk));
            
            % check too many bad entries are detected for each spectrum.
            badspc = (sum(logYif_isnan,1)/B) > th_badspc;
            logYif_isnan = or(logYif_isnan,badspc);
            
            lambda_r = 1./mad_expected.*(Ymdl)./(B*20);
            % do we need to set lambda_r to zero?? not sure.
            lambda_r(logYif_isnan) = 0;
    end

    lambda_c(logYif_isnan) = inf; lambda_c(~logYif_isnan) = 0;
    lambda_c([1,Nc],:,:) = 0; % safeguard
    
    
    resNrm = nansum(abs(lambda_r .* RR),[1,2]);
    
    % update logt_est
    if ffc_mode
        RR = RR + logT;
        Xtc = ones(1,L,S,precision,gpu_varargin);
    else
        if batch
            RR  = RR + pagefun(@mtimes,A(:,1,:),X(1,:,:));
        else
            RR  = RR + A(:,1)*X(1,:);
        end
        %Xtc = X(idxAlogT,:,:);
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
        % RR = RR - pagefun(@mtimes,logt_est,Xtc);
    else
        RR = RR - logt_est*X(1,:);
        % RR = RR - logt_est*Xtc;
    end
    resNewNrm = nansum(abs(lambda_r .* RR),[1,2]);
    
    % do we need to update lambda here? not sure.
    
    A(:,idxAlogT,:) = logt_est;
    
    switch upper(lambda_update_rule)
        case 'L1SUM'
            lambda_a_2(2:end,:,:) = lambda_a_2(2:end,:,:) .* resNewNrm ./ resNrm;
        case 'MED'
            error('not implemented yet');
        case 'NONE'
        otherwise
            error('Undefined LAMBDA_UPDATE_RULE: %s',lambda_update_rule);
    end
    
    
end

%%
%-------------------------------------------------------------------------%
% last iteration
%-------------------------------------------------------------------------%
if ffc_mode
    if batch
        [ X,Z,D,rho,Rhov,~,~,cost_val ] = huwacbl1_admm_gat_a_batch(A,...
            logYif-logt_est,C,...
            'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
            'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,...%'rho',rho,'Rhov',Rhov,...
            'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
            'PRECISION',precision,'YNORMALIZE',y_normalize);
    else
         [ X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
            = huwacbl1_admm_gat_a(A,logYif-logt_est,WA,'LAMBDA_A',lambda_a_2,...
            'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
            'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
            'tol',tol_huwacb,'maxiter',maxiter_huwacb,'verbose','no',...
            'precision',precision,'gpu',gpu,'Concavebase',C,...
            'debug',debug_huwacb,'YNORMALIZE',y_normalize);
    end
else
    if batch
        [ X,Z,D,rho,Rhov,~,~,cost_val ] = huwacbl1_admm_gat_a_batch(A,...
            logYif,C,...
            'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
            'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,...%'rho',rho,'Rhov',Rhov,...
            'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
            'PRECISION',precision,'YNORMALIZE',y_normalize);
    else
        [ X,Z,~,~,D,rho,Rhov,~,~,cost_val]...
            = huwacbl1_admm_gat_a(A,logYif,WA,'LAMBDA_A',lambda_a_2,...
            'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
            'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
            'tol',tol_huwacb,'maxiter',maxiter_huwacb,'verbose','no',...
            'precision',precision,'gpu',gpu,'Concavebase',C,...
            'debug',debug_huwacb,'YNORMALIZE',y_normalize);
    end
end

if ffc_mode     
    Xtc = ones(1,L,S,precision,gpu_varargin);
else
    Xtc = X(idxAlogT,:,:);
end

% substituting all the variables
if batch
    logAB        = pagefun(@mtimes,Alib,X(idxAlib,:,:));
    logBg        = pagefun(@mtimes,C,Z);
    logIce       = pagefun(@mtimes,Aicelib,X(idxAice,:,:));
    logYif_cor   = logYif - pagefun(@mtimes,logt_est,Xtc) - logIce;
else
    logAB = Alib*X(idxAlib,:);
    logBg = C*Z;
    logIce = Aicelib*X(idxAice,:);
    logYif_cor = logYif - logt_est*Xtc - logIce;
end
% logt_est     = A(:,idx,:);
logYif_cor(logYif_isnan) = nan;

if batch
    [logYif_cor,logt_est,logAB,logBg,logIce,logYif_isnan,X,badspc]...
        = gather(logYif_cor,logt_est,logAB,logBg,logIce,logYif_isnan,X,badspc);
else
end

Xt           = Xtc;
Xice         = X(idxAice,:,:);
Xlib         = X(idxAlib,:,:);

end

