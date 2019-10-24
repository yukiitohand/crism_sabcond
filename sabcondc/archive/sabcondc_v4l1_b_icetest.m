% Aice will be kept until the end
function [ water_ice_exist, Xice1_mean_insig, Xice1, idx_insig,ancillary] = sabcondc_v4l1_b_icetest( Alib,Aice,logYifc,wvc,logtc,varargin )
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
maxiter_lad = 1000;
tol_lad = 1e-4;
verbose_lad = 'no';
nIter = 5;
vis = 0;
lambda_a = 0.01;
lam_a_ice = 0;
logYifc_cat = 0;
logYraifc_cat = 0;
verbose_huwacb = 'no';
isdebug = false;
mode_name = 'SingleStep'; % {'SingleStep','DoubleStep'}
spike_detection_opt = 'original'; % {'original' or 'incremental'}
gp = [];
threshold_insig_logAB_bprmvd = 0.5;
threshold_Xice1 = 0.1;

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
            case 'LAMBDA_A_ICE'
                lam_a_ice = varargin{i+1};
            case 'LOGYIFC_CAT'
                logYifc_cat = varargin{i+1};
            case 'LOGYRAIFC_CAT'
                logYraifc_cat = varargin{i+1};
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

nIce = size(Aice,2); nlib = size(Alib,2); nlogtc = size(logtc,2);

A = [logtc Aice Alib];
nA1 = size(A,2);

idxAicestrt = nlogtc+1;
idxAlibstrt = nlogtc+nIce+1;
idxAlogtc = false(1,size(A,2));
idxAlogtc(1:idxAicestrt-1) = true;
idxAice = false(1,size(A,2));
idxAice(idxAicestrt:idxAlibstrt-1) = true;
idxAlib = false(1,size(A,2));
idxAlib(idxAlibstrt:end) = true;


gp_bool = (gp==1);
% bp_nan = isnan(gp);
logYifc_bprmvd = logYifc(gp_bool,:);
logYifc_bprmvd_ori = logYifc_bprmvd; % save for later
wvc_bprmvd = wvc(gp_bool,:);
A_bprmvd = A(gp_bool,:);
Alib_bprmvd = Alib(gp_bool,:);
if ~isempty(Aice)
    Aice_bprmvd = Aice(gp_bool,:);
end
% [ groups ] = nan_grouper( logYifc_ori );
logYifc_bprmvd_isnan = isnan(logYifc_bprmvd);
logYifc_bprmvd_isnan_ori = logYifc_bprmvd_isnan; % save for later
% in case there is any nan (corresponding to negative values in the original domain)
% no model is learned, so just interpolate with neighboring spectels
logYifc_bprmvd = interp_nan_column(logYifc_bprmvd,logYifc_bprmvd_isnan,wvc_bprmvd);
[L_bprmvd,Ny] = size(logYifc_bprmvd);

vldpxl = (sum(logYifc_bprmvd_isnan,1)/L_bprmvd) < 0.8;

%% peform unmixing
lambda_a_1 = zeros([1,nA1]);
lambda_tmp = lambda_a;
if ~isempty(Aice)
    lam_tmp_ice = lam_a_ice;
    lambda_a_1(idxAice) = lam_tmp_ice;
end
lambda_a_1(idxAlib) = lambda_tmp;

X1 = nan(nA1,Ny);
Z1 = nan(L_bprmvd,Ny);
D1 = nan([nA1+L_bprmvd*2,Ny]);

[ X1(:,vldpxl),Z1(:,vldpxl),C1,~,D1(:,vldpxl),rho,Rhov,~]...
    = huwacbl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,'LAMBDA_A',lambda_a_1,...
                               'tol',1e-5,'maxiter',maxiter_huwacb,'verbose',verbose_huwacb);
logBg_bprmvd = C1*Z1;
logYifc_model_bprmvd = logBg_bprmvd + A_bprmvd*X1;


switch spike_detection_opt
%     case 'incremental'
%         rr1 = logYifc_bprmvd-logYifc_model_bprmvd;
%         logYifc_bprmvd_isnan = or(logYifc_bprmvd_isnan,abs(rr1)>0.05);
%         logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
    case 'original'
        rr1 = logYifc_bprmvd_ori-logYifc_model_bprmvd;
        rr1_std = std(rr1,[],2);
        logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(rr1)>0.1 );
        logYifc_bprmvd_isnan(rr1_std<0.015,:) = false; % if the standard deviation of the channel is small enough, bring back to good one.
        logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
        vldpxl = ~any(isnan(logYifc_bprmvd),1);
%     case 'so_original'
%         rr1 = logYifc - logYifc_model_bprmvd;
%         logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(rr1)>0.05 );
    otherwise 
        error('not implemented');
end

if isdebug
    %figure(10); imsc(logYifc_bprmvd_isnan); 
    %title('Bad pixels'); 
    %drawnow;
end

Xlib = X1(idxAlib,:); Xlogtc = X1(idxAlogtc,:);
Xice = X1(idxAice,:);
logAB_bprmvd = Alib_bprmvd*Xlib;
logIce_bprmvd = Aice_bprmvd*Xice;

%% test ice
% select bland spectra
idx_insig = abs(sum(logAB_bprmvd,1))<threshold_insig_logAB_bprmvd;

% amount of ice
Xice1 = sum(Xice,1);
Xice1_mean_insig = mean(Xice1(idx_insig));

water_ice_exist = Xice1_mean_insig > threshold_Xice1;


ancillary = [];
ancillary.X = X1;
ancillary.lambda = [];
ancillary.lambda.init = lambda_a; % edited by Yuki 2017/10/11
%ancillary.lambda.last = lambda_a_2(2:end); % edited by Yuki 2017/10/11, the final lambda used
ancillary.nIter = nIter;
ancillary.huwacb_func = 'huwacbl1_gadmm_a_v2';
ancillary.maxiter_huwacb = maxiter_huwacb;
ancillary.tol_huwacb = tol_huwacb;
ancillary.gp_bool = gp_bool;


end