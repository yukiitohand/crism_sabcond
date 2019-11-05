function [ logt_est,logYifc_cor,logAB,logBg,logIce,logYifc_cor_ori,logYifc_isnan,...
    ancillary,vldpxl]...
    = sabcondc_v3l1_pub( Alib,logYifc,wvc,logtc,varargin )
% [ logt_est,logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,...
%    ancillary,vldpxl]...
%    = sabcondc_v3l1_pub( Alib,logYifc,wvc,logtc,varargin )
% Perfrom atmospheric and de-noising of CRISM data with Algorithm 1 (bad
% entries are replaced with model values at each iteration).
%
% INPUT PARAMETERS
%   Alib: array, [B x Nlib]
%        library matrix
%   logYifc: array, [B x L]
%       input I/F spectra of one spatial column
%   wvc: array, [B x 1]
%       wavelength frame
%   logtc: array, [B x Ntc]
%       collection of transmission spectra
%
% OUTPUT PARAMETERS
%   logt_est: array, [B x 1]
%       estimated transmission spectrum
%   logYifc_cor: array, [B x L]
%       corrected i/f spectra
%   logAB: array, [B x L]
%       estimated absorption spectra
%   logBg: array, [B x L]
%       estimated background spectra
%   logIce: array, [B x L]
%       estimated ice contributions
%   logYifc_cor_ori: array, [B x L]
%       corrected i/f spectra (bad entries are also processed)
%   logYifc_isnan: boolean array, [B x L]
%       bad entries are flagged.
%   ancillary: ancillary information
%       Field contains
%           X: [(1+Nlib+Nice) x L] estimated abundance matrix
%           gp_bool: boolean array, [L x 1]
%   vldpxl: boolean array, [1 x L]
%       flag if the spectra has sufficiently a small number of bad entries
%       (<THREHOLD_BADSPC)
%
% OPTIONAL PARAMETERS
%  ## GENERAL PARAMETERS #-------------------------------------------------
%   'AICELIB': array, [L x Nc]
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
%        (default) false
%   'PRECISION': {'single','doulbe'}
%        data type used for calculation
%        (default) 'double'

%% VARARGIN
% ## GENERAL PARAMETERS #--------------------------------------------------
Aicelib = [];
nIter = 5;
th_badspc = 0.8;
gp = [];
% ## HUWACB PARAMETERS #---------------------------------------------------
lambda_a = 0.01;
lambda_a_ice = 0;
maxiter_huwacb = 100;
tol_huwacb = 1e-4;
verbose_huwacb = 'no';
debug_huwacb = false;
% ## LAD PARAMETERS #------------------------------------------------------
maxiter_lad = 1000;
tol_lad = 1e-4;
verbose_lad = 'no';
debug_lad = false;
% ## PROCESSING PARAMETERS #-----------------------------------------------
gpu = false;
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
            case 'GP'
                gp = varargin{i+1};
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
            case 'DEBUG_HUWACB'
                debug_huwacb = varargin{i+1};
            case 'VERBOSE_HUWACB'
                verbose_huwacb = varargin{i+1};
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
            case 'PRECISION'
                precision = varargin{i+1};
            
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

[B,Ny]   = size(logYifc);
[~,Nlib] = size(Alib);
[~,Ntc]  = size(logtc);
[~,Nice] = size(Aicelib);

if Nice==0
    Aicelib = zeros(B,0,precision);
end

gp_bool = (gp==1);
B_bprmvd = sum(gp_bool);
logYifc_bprmvd = logYifc(gp_bool,:);
logYifc_bprmvd_ori = logYifc_bprmvd; % save for later
wvc_bprmvd = wvc(gp_bool,:);
Alib_bprmvd = Alib(gp_bool,:);
logtc_bprmvd = logtc(gp_bool,:);
Aice_bprmvd = Aicelib(gp_bool,:);

A_bprmvd = [logtc_bprmvd Aice_bprmvd Alib_bprmvd];
N_A1 = Ntc + Nice + Nlib;

idxAicestrt = Ntc+1;
idxAlibstrt = Ntc+Nice+1;
idxAlogtc = false(1,N_A1);
idxAlogtc(1:idxAicestrt-1) = true;
idxAice = false(1,N_A1);
idxAice(idxAicestrt:idxAlibstrt-1) = true;
idxAlib = false(1,N_A1);
idxAlib(idxAlibstrt:end) = true;

logYifc_bprmvd_isnan = isnan(logYifc_bprmvd);
logYifc_bprmvd_isnan_ori = logYifc_bprmvd_isnan; % save for later
% in case there is any nan (corresponding to negative values in the original domain)
% no model is learned, so just interpolate with neighboring spectels
logYifc_bprmvd = interp_nan_column(logYifc_bprmvd,logYifc_bprmvd_isnan,wvc_bprmvd);
vldpxl = (sum(logYifc_bprmvd_isnan,1)/B_bprmvd) < 0.5;

%% initialization of logt_est
lambda_a_1 = zeros(N_A1,1,precision);
lambda_tmp = lambda_a;
lambda_tmp_ice = lambda_a_ice;
% lambda_a_1(idxAlibstrt:end) = lambda_tmp;
lambda_a_1(idxAlib) = lambda_tmp;
lambda_a_1(idxAice) = lambda_tmp_ice;

X1 = nan(N_A1,Ny,precision); Z1 = nan(B_bprmvd,Ny,precision); 
D1 = nan(N_A1+B_bprmvd*2,Ny,precision);

[ X1(:,vldpxl),Z1(:,vldpxl),C1,~,D1(:,vldpxl),rho,Rhov,~,~,cost_val] ...
    = huwacbl1_admm_gat_a(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
    'LAMBDA_A',lambda_a_1,'tol',1e-5,'maxiter',maxiter_huwacb,...
    'verbose',verbose_huwacb,'precision',precision,'gpu',gpu,'debug',debug_huwacb);

logYifc_model_bprmvd = A_bprmvd*X1+C1*Z1;
% noise detection and replacement
RR_bprmvd = logYifc_bprmvd_ori-logYifc_model_bprmvd;
RR_std = nanstd(RR_bprmvd,[],2);
logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(RR_bprmvd)>0.1 );
logYifc_bprmvd_isnan(RR_std<0.015,:) = false; % if the standard deviation of the channel is small enough, bring back to good one.
logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
vldpxl = ~any(isnan(logYifc_bprmvd),1);

% Xlib = X1(idxAlib,:); Xlogtc = X1(idxAlogtc,:);
RR_bprmvd = logYifc_bprmvd - logYifc_model_bprmvd;
resNrm_bprmvd = nansum(abs(RR_bprmvd),'all');
RR_bprmvd = RR_bprmvd + A_bprmvd(:,idxAlogtc) * X1(idxAlogtc,:);

Xlogtc_1d = sum(X1(idxAlogtc,:),1);
% r_lad = zeros(Ny,B_bprmvd,precision); d_lad = zeros(Ny+1,B_bprmvd,precision);
Rhov_lad = ones(1+Ny,1,precision);
[logt_est_bprmvd,~,~,rho_lad,Rhov_lad([true vldpxl],:),~,~,cost_val]...
    = lad_admm_gat_b(Xlogtc_1d(:,vldpxl)', RR_bprmvd(:,vldpxl)',...
             'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad,...
             'PRECISION',precision,'gpu',gpu,'debug',debug_lad);

logt_est_bprmvd = logt_est_bprmvd';
RR_bprmvd = RR_bprmvd - logt_est_bprmvd*Xlogtc_1d;
resNewNrm_bprmvd = nansum(abs(RR_bprmvd),'all');

%%
%-------------------------------------------------------------------------%
% main loop
%-------------------------------------------------------------------------%
A_bprmvd = [logt_est_bprmvd Aice_bprmvd Alib_bprmvd];
N_A = size(A_bprmvd,2);
X = [Xlogtc_1d;X1(or(idxAice,idxAlib),:)];
% D = [sum(D1(idxAlogtc,:),1);D1(idxAlib,:)];
D = [zeros(1,size(D1,2),precision); D1(idxAicestrt:end,:)];
C = C1;
Z = Z1;


rho = ones([1,Ny]);

clear X1 Z1 C1 D1;

idxAice = false(1,N_A); idxAice(2:(Nice+1)) = true;
idxAlib = false(1,N_A); idxAlib((Nice+2):end) = true;

% always update lambda_tmp
lambda_a_2 = zeros(N_A,1);
lambda_tmp = lambda_tmp*resNewNrm_bprmvd/resNrm_bprmvd;
lambda_tmp_ice = lambda_tmp_ice*resNewNrm_bprmvd/resNrm_bprmvd;

for j=2:nIter+1
    % lambda_a_2(2:end) = lambda_tmp;   
    % rr = logYifc_bprmvd - A_bprmvd*X - C*Z;
    lambda_a_2(idxAlib) = lambda_tmp;
    lambda_a_2(idxAice) = lambda_tmp_ice;
    if j==2
        [ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
            = huwacbl1_admm_gat_a(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
                            'LAMBDA_A',lambda_a_2,'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
                            'D0',D(:,vldpxl),'X0',X(:,vldpxl),...
                            'R0',RR_bprmvd(:,vldpxl),...
                            'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb,...
                            'PRECISION',precision,'gpu',gpu,'debug',debug_huwacb);
    else
       [ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
           = huwacbl1_admm_gat_a(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
                            'LAMBDA_A',lambda_a_2,'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
                            'D0',D(:,vldpxl),'X0',X(:,vldpxl),...
                            'R0',RR_bprmvd(:,vldpxl),'rho',rho(:,vldpxl),'Rhov',Rhov,...
                            'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
                            'PRECISION',precision,'gpu',gpu,'debug',debug_huwacb);
    end
    
    %logBg_bprmvd = C*Z;
    %logAB_bprmvd = Alib_bprmvd*X(2:end,:);
    logYifc_model_bprmvd = A_bprmvd*X + C*Z;
    
    % spike detection
    RR_bprmvd = logYifc_bprmvd_ori-logYifc_model_bprmvd;
    logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(RR_bprmvd)>0.015 );
    logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
    vldpxl = (sum(logYifc_bprmvd_isnan,1)/B_bprmvd) < th_badspc;
    
    %evaluate residual again
    RR_bprmvd = logYifc_bprmvd - logYifc_model_bprmvd;
    resNrm_bprmvd = nansum(abs(RR_bprmvd),'all');

    %update logt_est!
    RR_bprmvd = RR_bprmvd + A_bprmvd(:,1)*X(1,:);
    [logt_est_bprmvd,~,~,...
        rho_lad,Rhov_lad([true vldpxl],:),~,~,cost_val]...
        = lad_admm_gat_b(X(1,vldpxl)', RR_bprmvd(:,vldpxl)',...
        'rho',rho_lad,'Rhov',Rhov_lad([true vldpxl],:),...
        'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad,...
        'PRECISION',precision,'gpu',gpu,'debug',debug_lad);
    
    logt_est_bprmvd = logt_est_bprmvd';
    
    %
    RR_bprmvd = RR_bprmvd - logt_est_bprmvd*X(1,:);
    resNewNrm_bprmvd = nansum(abs(RR_bprmvd),'all');
    
    A_bprmvd(:,1) = logt_est_bprmvd;
    
    lambda_tmp = lambda_tmp*resNewNrm_bprmvd/resNrm_bprmvd;
    lambda_tmp_ice = lambda_tmp_ice*resNewNrm_bprmvd/resNrm_bprmvd;

end

%%
%-------------------------------------------------------------------------%
%last iteration after estimating log_est
%-------------------------------------------------------------------------%
% rr = logYifc_bprmvd - A_bprmvd*X - C*Z;
[ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
    = huwacbl1_admm_gat_a(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
                            'LAMBDA_A',lambda_a_2,'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
                            'D0',D(:,vldpxl),'X0',X(:,vldpxl),'R0',RR_bprmvd(:,vldpxl),...
                            'rho',rho(:,vldpxl),'Rhov',Rhov,...
                            'verbose',verbose_huwacb,...
                            'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
                            'PRECISION',precision,'gpu',gpu,'debug',debug_huwacb);

%%
% original wavelength channels
logBg = nan(B,Ny,precision); logBg(gp_bool,:) = C*Z;
logAB = Alib*X(idxAlib,:); 
logIce = Aicelib*X(idxAice,:);

logt_est = nan(B,1,precision); logt_est(gp_bool) = A_bprmvd(:,1);

% corrected spectra
logYifc_cor = nan(B,Ny,precision);
logYifc_bprmvd(logYifc_bprmvd_isnan) = nan;
logYifc_cor(gp_bool,:) = logYifc_bprmvd - A_bprmvd(:,1)*X(1,:)-logIce(gp_bool,:);

% corrected spectra
logYifc_cor_ori = nan(B,Ny,precision);
logYifc_cor_ori(gp_bool,:) = logYifc_bprmvd_ori - A_bprmvd(:,1)*X(1,:) - logIce(gp_bool);

logYifc_isnan = true(B,Ny);
logYifc_isnan(gp_bool,:) = logYifc_bprmvd_isnan;

% 
logBg = interp_nan_column(logBg,logYifc_isnan,wvc);

vldpxl = (sum(logYifc_bprmvd_isnan,1)/B_bprmvd) < th_badspc;

ancillary = [];
ancillary.X = X;
ancillary.gp_bool = gp_bool;


end