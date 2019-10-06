function [ logt_est,logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,...
    ancillary,rr_ori,vldpxl]...
    = sabcondc_v3l1_pub( Alib,logYifc,wvc,logtc,varargin )
% [ logt_est,logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,...
%    ancillary,rr_ori,vldpxl]...
%    = sabcondc_v3l1_pub( Alib,logYifc,wvc,logtc,varargin )
% Perfrom atmospheric and de-noising of CRISM data using 

%%
%-------------------------------------------------------------------------%
% Variable setting
%-------------------------------------------------------------------------%
maxiter_huwacb = 100;
% tol_huwacb = 1e-4;
tol_huwacb = 1e-4;
maxiter_lad = 1000;
tol_lad = 1e-4;
verbose_lad = 'no';
debug_lad = false;
nIter = 5;
lambda_a = 0.01;
verbose_huwacb = 'no';
debug_huwacb = false;
gp = [];
precision = 'double';
gpu = false;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'GP'
                gp = varargin{i+1};
            case 'NITER'
                nIter = varargin{i+1};
            case 'LAMBDA_A'
                lambda_a = varargin{i+1};
            case 'MAXITER_HUWACB'
                maxiter_huwacb = round(varargin{i+1});
                if (maxiter_huwacb <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'TOL_HUWACB'
                tol_huwacb = varargin{i+1};
            case 'VERBOSE_HUWACB'
                verbose_huwacb = varargin{i+1};
            case 'MAXITER_LAD'
                maxiter_lad = round(varargin{i+1});
                if (maxiter_huwacb <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'TOL_LAD'
                tol_lad = varargin{i+1};
            case 'VERBOSE_LAD'
                verbose_lad = varargin{i+1};
            case 'DEBUG_HUWACB'
                debug_huwacb = varargin{i+1};
            case 'DEBUG_LAD'
                debug_lad = varargin{i+1};
            case 'PRECISION'
                precision = varargin{i+1};
            case 'GPU'
                gpu = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

[B,Ny] = size(logYifc);
[~,Nlib] = size(Alib);
[~,Ntc]  = size(logtc);


gp_bool = (gp==1);
B_bprmvd = sum(gp_bool);
logYifc_bprmvd = logYifc(gp_bool,:);
logYifc_bprmvd_ori = logYifc_bprmvd; % save for later
wvc_bprmvd = wvc(gp_bool,:);
Alib_bprmvd = Alib(gp_bool,:);
logtc_bprmvd = logtc(gp_bool,:);
A_bprmvd = [logtc_bprmvd Alib_bprmvd];
N_A1 = Ntc + Nlib;

idxAlibstrt = Ntc+1;
idxAlogtc = false(1,N_A1);
idxAlogtc(1:idxAlibstrt-1) = true;
idxAlib = false(1,N_A1);
idxAlib(idxAlibstrt:end) = true;

logYifc_bprmvd_isnan = isnan(logYifc_bprmvd);
logYifc_bprmvd_isnan_ori = logYifc_bprmvd_isnan; % save for later
% in case there is any nan (corresponding to negative values in the original domain)
% no model is learned, so just interpolate with neighboring spectels
logYifc_bprmvd = interp_nan_column(logYifc_bprmvd,logYifc_bprmvd_isnan,wvc_bprmvd);
vldpxl = (sum(logYifc_bprmvd_isnan,1)/B_bprmvd) < 0.8;

%% initialization of logt_est
lambda_a_1 = zeros(1,N_A1,precision);
lambda_tmp = lambda_a;
% lambda_a_1(idxAlibstrt:end) = lambda_tmp;
lambda_a_1(idxAlib) = lambda_tmp;

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
A_bprmvd = [logt_est_bprmvd Alib_bprmvd];
X = [Xlogtc_1d;X1(idxAlib,:)];
% D = [sum(D1(idxAlogtc,:),1);D1(idxAlib,:)];
D = [zeros(1,size(D1,2),precision); D1(idxAlibstrt:end,:)];
C = C1;
Z = Z1;
lambda_a_2 = zeros(1,1+Nlib);

rho = ones([1,Ny]);

clear X1 Z1 C1 D1;

% always update lambda_tmp
lambda_tmp = lambda_tmp*resNewNrm_bprmvd/resNrm_bprmvd;

for j=2:nIter+1
    lambda_a_2(2:end) = lambda_tmp;   
    % rr = logYifc_bprmvd - A_bprmvd*X - C*Z;
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
    vldpxl = (sum(logYifc_bprmvd_isnan,1)/B_bprmvd) < 0.8;
    
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
logAB = Alib*X(2:end,:);

logt_est = nan(B,1,precision); logt_est(gp_bool) = A_bprmvd(:,1);

% corrected spectra
logYifc_cor = nan(B,Ny,precision);
logYifc_bprmvd(logYifc_bprmvd_isnan) = nan;
logYifc_cor(gp_bool,:) = logYifc_bprmvd - A_bprmvd(:,1)*X(1,:);

% corrected spectra
logYifc_cor_ori = nan(B,Ny,precision);
logYifc_cor_ori(gp_bool,:) = logYifc_bprmvd_ori - A_bprmvd(:,1)*X(1,:);

logYifc_isnan = true(B,Ny);
logYifc_isnan(gp_bool,:) = logYifc_bprmvd_isnan;

% 
logBg = interp_nan_column(logBg,logYifc_isnan,wvc);

% residual
rr_ori = nan(B,Ny,precision);
rr_ori(gp_bool,:) = logYifc_bprmvd_ori - logAB(gp_bool,:) - logBg(gp_bool,:) - A_bprmvd(:,1)*X(1,:);

vldpxl = (sum(logYifc_bprmvd_isnan,1)/B_bprmvd) < 0.7;


ancillary = [];
ancillary.X = X;
ancillary.lambda = [];
ancillary.lambda.init = lambda_a; % edited by Yuki 2017/10/11
ancillary.lambda.last = lambda_a_2(2:end); % edited by Yuki 2017/10/11, the final lambda used
ancillary.nIter = nIter;
ancillary.huwacb_func = 'huwacbl1_gadmm_a_v2';
ancillary.maxiter_huwacb = maxiter_huwacb;
ancillary.tol_huwacb = tol_huwacb;
ancillary.gp_bool = gp_bool;


end