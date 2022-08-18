function [ logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,ancillary,rr_ori,vldpxl] = sabcondc_v3l1_givent( Alib,logYifc,wvc,logt_est,varargin )
% [ logt_est ] = init_sabcond( Alib,logYifc,wvc,logtc,t_mode )
%   initialize the transmission spectrum using the transmission spectra of
%   ADR data record
%   Input Parameters
%     Alib: library matrix [L x Na]
%     logYifc: log reflectance [L x N]
%     wvc: wavelength [L x 1]
%     logtc: log of transmission spectrum in the ADR [L x ..], the number
%            of columns is 1 when t_mode=1 and multiple when t_mode = 2
%     t_mode: {1,2},1: one transmission, 2: collection of transmission (removed obsoleted)
%   Optional Parameters
%     REFINEMENT: integer, number of refinement iterations
%                 (default) 0
%   Output Parameters
%     logt_est: initialized transmission spectrum [L x 1]

maxiter_huwacb = 100;
% tol_huwacb = 1e-4;
tol_huwacb = 1e-4;
nIter = 5;
vis = 0;
lambda_a = 0.01;
verbose_huwacb = 'no';
isdebug = false;
mode_name = 'SingleStep'; % {'SingleStep','DoubleStep'}
spike_detection_opt = 'original'; % {'original' or 'incremental'}
gp = [];

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'GP'
                gp = varargin{i+1};
            case 'NITER'
                nIter = varargin{i+1};
            case 'VIS'
                vis = varargin{i+1};
            case 'LAMBDA_A'
                lambda_a = varargin{i+1};
            case 'MAXITER'
                maxiter_huwacb = round(varargin{i+1});
                if (maxiter_huwacb <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'TOL'
                tol_huwacb = varargin{i+1};
            case 'MODE'
                mode_name = varargin{i+1};
            case 'SPIKE_DETECTION_OPT'
                spike_detection_opt = varargin{i+1};
            case 'VERBOSE'
                verbose_huwacb = varargin{i+1};
            case 'DEBUG'
                isdebug = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

A = [logt_est Alib];
NA = size(A,2);

gp_bool = (gp==1);
bp_nan = isnan(gp);
logYifc_bprmvd = logYifc(gp_bool,:);
logYifc_bprmvd_ori = logYifc_bprmvd; % save for later
wvc_bprmvd = wvc(gp_bool,:);
A_bprmvd = A(gp_bool,:);
Alib_bprmvd = Alib(gp_bool,:);
% [ groups ] = nan_grouper( logYifc_ori );
logYifc_bprmvd_isnan = isnan(logYifc_bprmvd);
logYifc_bprmvd_isnan_ori = logYifc_bprmvd_isnan; % save for later
% in case there is any nan (corresponding to negative values in the original domain)
% no model is learned, so just interpolate with neighboring spectels
logYifc_bprmvd = interp_nan_column(logYifc_bprmvd,logYifc_bprmvd_isnan,wvc_bprmvd);

[L_bprmvd,Ny] = size(logYifc_bprmvd);
vldpxl = (sum(logYifc_bprmvd_isnan,1)/L_bprmvd) < 0.8;

lambda_tmp = lambda_a;



%% perform un

X = nan(NA,Ny);
% D = [sum(D1(idxAlogtc,:),1);D1(idxAlib,:)];
D = nan(NA+L_bprmvd*2,Ny);
Z = nan(L_bprmvd,L_bprmvd);
lambda_a_2 = zeros([1,size(A_bprmvd,2)]);
rho = ones([1,Ny]);

% always update lambda_tmp

lambda_a_2(2:end) = lambda_tmp;  

[ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
            = huwacbl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
                            'LAMBDA_A',lambda_a_2,...
                            'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb);

logBg_bprmvd = C*Z;
logAB_bprmvd = Alib_bprmvd*X(2:end,:);
logYifc_model_bprmvd = logBg_bprmvd + logAB_bprmvd + A_bprmvd(:,1)*X(1,:);
    
switch spike_detection_opt
    case 'original'
        rr = logYifc_bprmvd_ori-logYifc_model_bprmvd;
        logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(rr)>0.015 );
        logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
        vldpxl = (sum(logYifc_bprmvd_isnan,1)/L_bprmvd) < 0.8;
    otherwise 
        error('not implemented');
end

R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
res_bprmvd = R_bprmvd-A_bprmvd(:,1)*X(1,:);
resNrm_bprmvd = nansum(nansum(abs(res_bprmvd)));
resNew_bprmvd = R_bprmvd - A_bprmvd(:,1)*X(1,:);
resNewNrm_bprmvd = nansum(nansum(abs(resNew_bprmvd)));
    
lambda_tmp = lambda_tmp*resNewNrm_bprmvd/resNrm_bprmvd;




%% last iteration after estimating log_est
rr = logYifc_bprmvd - A_bprmvd*X - C*Z;
[ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
    = huwacbl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
                            'LAMBDA_A',lambda_a_2,'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
                            'D0',D(:,vldpxl),'X0',X(:,vldpxl),'R0',rr(:,vldpxl),...
                            'rho',rho(:,vldpxl),'Rhov',Rhov,...
                            'verbose',verbose_huwacb,...
                            'tol',tol_huwacb,'maxiter',maxiter_huwacb);

% original wavelength channels
logBg = nan(size(logYifc));
logBg(gp_bool,:) = C*Z;

logAB = Alib*X(2:end,:);

logt_est = nan(size(logYifc,1),1);
logt_est(gp_bool) = A_bprmvd(:,1);

% corrected spectra
logYifc_cor = nan(size(logYifc));
logYifc_bprmvd(logYifc_bprmvd_isnan) = nan;
logYifc_cor(gp_bool,:) = logYifc_bprmvd - A_bprmvd(:,1)*X(1,:);
% corrected spectra
logYifc_cor_ori = nan(size(logYifc));
logYifc_cor_ori(gp_bool,:) = logYifc_bprmvd_ori - A_bprmvd(:,1)*X(1,:);

logYifc_isnan = true(size(logYifc));
logYifc_isnan(gp_bool,:) = logYifc_bprmvd_isnan;

logBg = interp_nan_column(logBg,logYifc_isnan,wvc);

% residual
rr_ori = nan(size(logYifc));
rr_ori(gp_bool,:) = logYifc_bprmvd_ori - logAB(gp_bool,:) - logBg(gp_bool,:) - A_bprmvd(:,1)*X(1,:);

vldpxl = (sum(logYifc_bprmvd_isnan,1)/L_bprmvd) < 0.7;


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

