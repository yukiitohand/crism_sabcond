function [ logt_est,logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,ancillary,rr_ori,vldpxl,res_exp,std_expected_scaled] = sabcondc_v5l1_med( Alib,logYifc,wvc,logtc,varargin )
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

maxiter_huwacb = 200;
% tol_huwacb = 1e-4;
tol_huwacb = 1e-4;
verbose_huwacb = false;
maxiter_lad = 1000;
tol_lad = 1e-5;
verbose_lad = false;
nIter = 5;
vis = 0;
lambda_a = 0.01;
logYifc_cat = [];
logYraifc_cat = [];
stdl1_ifdf = '';
Yif = [];
T = [];
BPpri1nan = [];
BPall1nan = [];
photon_mad = [];
lbl = [];
SFimgc = [];
WA_um_pitch = [];

isdebug = false;
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
            case 'STDL1_IFDF'
                stdl1_ifdf = varargin{i+1};
            case 'YIF'
                Yif = varargin{i+1};
            case 'T'
                T = varargin{i+1};
            case 'TOL'
                tol_huwacb = varargin{i+1};
            case 'SPIKE_DETECTION_OPT'
                spike_detection_opt = varargin{i+1};
            case 'VERBOSE'
                verbose_huwacb = varargin{i+1};
            case 'DEBUG'
                isdebug = varargin{i+1};
            case 'LOGYIFC_CAT'
                logYifc_cat = varargin{i+1};
            case 'LOGYRAIFC_CAT'
                logYraifc_cat = varargin{i+1};
            case 'MAXITER'
                maxiter_huwacb = round(varargin{i+1});
                if (maxiter_huwacb <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'BP_PRI'
                BPpri1nan = varargin{i+1};
            case 'BP_ALL'
                BPall1nan = varargin{i+1};
            case 'PHOTON_MAD'
                photon_mad = varargin{i+1};
            case 'LBL'
                lbl = varargin{i+1};
            case 'SFIMGC'
                SFimgc = varargin{i+1};
            case 'WA_UM_PITCH'
                WA_um_pitch = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

A = [logtc Alib];
idxAlibstrt = size(logtc,2)+1;
idxAlogtc = false(1,size(A,2));
idxAlogtc(1:idxAlibstrt-1) = true;
idxAlib = false(1,size(A,2));
idxAlib(idxAlibstrt:end) = true;

N_A1 = size(A,2); Ny = size(logYifc,2);


gp_bool = (gp==1);
logYifc_bprmvd = logYifc(gp_bool,:);
logYifc_bprmvd_ori = logYifc_bprmvd; % save for later
wvc_bprmvd = wvc(gp_bool,:);
A_bprmvd = A(gp_bool,:);
Alib_bprmvd = Alib(gp_bool,:);
SFimgc_bprmvd = SFimgc(gp_bool,:);
WA_um_pitch_bprmvd = WA_um_pitch(gp_bool);
logYifc_bprmvd_isnan = isnan(logYifc_bprmvd);
logYifc_bprmvd_isnan_ori = logYifc_bprmvd_isnan; % save for later
% in case there is any nan (corresponding to negative values in the original domain)
% no model is learned, so just interpolate with neighboring spectels
logYifc_bprmvd = interp_nan_column(logYifc_bprmvd,logYifc_bprmvd_isnan,wvc_bprmvd);
Yifc_bprmvd = interp_nan_column(Yif(gp_bool,:),logYifc_bprmvd_isnan,wvc_bprmvd);

logYifc_bprmvd_good_1nan = double(~logYifc_bprmvd_isnan);
logYifc_bprmvd_good_1nan(logYifc_bprmvd_isnan) = nan;
logYifc_bprmvd_good_1nan_ori = logYifc_bprmvd_good_1nan;

[L_bprmvd,Ny] = size(logYifc_bprmvd);
vldpxl = (sum(logYifc_bprmvd_isnan,1)/L_bprmvd) < 0.7;

% logBg_bprmvd = nan(size(logYifc));
%% initialization of logt_est
lambda_a_1 = zeros([size(A,2),1]);
lambda_tmp = lambda_a;
lambda_a_1(idxAlib) = lambda_tmp;

X1 = zeros(N_A1,Ny);
Z1 = zeros(L_bprmvd,Ny);
D1 = zeros([N_A1+L_bprmvd*2,Ny]);

% [ X1(:,vldpxl),Z1(:,vldpxl),C1,~,D1(:,vldpxl),rho,Rhov,~] = huwacbl1_gadmm_a_v2(...
%     A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,'LAMBDA_A',lambda_a_1,...
%     'tol',1e-5,'maxiter',maxiter_huwacb,'verbose',verbose_huwacb);
%
mYif = nanmean(Yifc_bprmvd,2);
% mm = nanmean(mYif);
mYif1 = mYif'*Yifc_bprmvd/norm(mYif,2)^2;
T_bprmvd = T(gp_bool,:);
T_bprmvd_mean = nanmean(T_bprmvd,2);

%lambda_r_bprmvd = 1./stdl1_ifdf(gp_bool).*(mYif*mYif1)/400;
lambda_r_bprmvd = 1./(stdl1_ifdf(gp_bool)+photon_mad(gp_bool,:)).*(mYif*mYif1)/(L_bprmvd*20);
% lambda_r_bprmvd = 1./stdl1_ifdf(gp_bool).*(mYif*mYif1).*T_bprmvd_mean./(L_bprmvd*40);
lambda_r_bprmvd(logYifc_bprmvd_isnan) = 0;
% lambda_r_bprmvd = lambda_r(gp_bool,:);
% lambda_r_bprmvd = lambda_r_bprmvd;
%lambda_r_bprmvd = lambda_r_bprmvd/L_bprmvd;

lambda_c_bprmvd = 1e-8 * ones(L_bprmvd,Ny);
lambda_c_bprmvd(logYifc_bprmvd_isnan) = 0.01;
lambda_c_bprmvd([1,L_bprmvd],:) = 0; % no weight for edges
% edge robust??
lambda_c_bprmvd([2:4,(L_bprmvd-3):(L_bprmvd-1)],:) = 0.5;
% mad_expected_bprmvd = (stdl1_ifdf(gp_bool)+photon_mad(gp_bool,:))./(mYif*mYif1);
% mad_expected_bprmvd_bad = mad_expected_bprmvd>0.01;
% lambda_c_bprmvd(2:6,:) = mad_expected_bprmvd_bad(2:6,:)*0.01;
% lambda_c_bprmvd((L_bprmvd-3):(L_bprmvd-1),:) = mad_expected_bprmvd_bad((L_bprmvd-3):(L_bprmvd-1),:)*0.01;



[ X1(:,vldpxl),Z1(:,vldpxl),C1,~,D1(:,vldpxl),rho,Rhov,~] = huwacbwl1_gadmm_a_v2(...
    A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
    'LAMBDA_A',lambda_a_1,'Lambda_c',lambda_c_bprmvd(:,vldpxl),'lambda_r',lambda_r_bprmvd(:,vldpxl),...
    'tol',1e-5,'maxiter',maxiter_huwacb,'verbose',verbose_huwacb,'YNormalize',false);

% lambda_c_bprmvd = 1e-8 * ones(size(C1,2),sum(vldpxl));
                           
logBg_bprmvd = C1*Z1;
% logYifc_model_bprmvd = logBg_bprmvd + A_bprmvd*X1;
Xlib = X1(idxAlib,:); Xlogtc = X1(idxAlogtc,:);
logAB_bprmvd = Alib_bprmvd*Xlib;
R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
res_bprmvd = R_bprmvd - A_bprmvd(:,idxAlogtc) * Xlogtc;
%resNrm_bprmvd = nansum(nansum(lambda_r_bprmvd.*abs(res_bprmvd)));
resNrm_bprmvd = nanmedian(lambda_r_bprmvd.*abs(res_bprmvd).*logYifc_bprmvd_good_1nan,'all');


switch spike_detection_opt
    case 'original'
        Ymodel_bprmvd = exp(logAB_bprmvd+logBg_bprmvd+logtc(gp_bool,:)*Xlogtc);
        res_exp1 = Yifc_bprmvd - Ymodel_bprmvd;
        mad_bprmvd_rr1 = robust_v3('med_abs_dev_from_med',res_exp1,2,'NOutliers',10);
        mad_bprmvd_log1 = robust_v3('med_abs_dev_from_med',res_bprmvd,2,'NOutliers',10);
        
        % First step of de-noising is bad pixel detection.
        % Our bad pixel detection is based on the median absolute deviation
        % from the median. If that value is greater than a threshold, then
        % that will be excluded.
        % EXCEPTION
        %  If the residual in the log domain has only a small variance,
        %  then such pixels will be not considered as a bad pixel. The
        %  pixel is likely to be the error of a multiplicative components
        %  in the calibration.
        logYifc_bprmvd_isnan_bp = false(size(Ymodel_bprmvd));
        bp_est_bool = and(mad_bprmvd_rr1>0.0015,mad_bprmvd_log1>0.005);
        logYifc_bprmvd_isnan_bp(bp_est_bool,:) = true;
        logYifc_bprmvd_isnan = or(logYifc_bprmvd_isnan_bp,logYifc_bprmvd_isnan_ori);
        % for detection of temporal spikes are not evaluated to avoid bias
        % caused by the initial transmission spectrum or inaccurate BK or
        % BI.
        
        % if the standard deviation of dark frames is small enough, bring back to good one.
%         logYifc_bprmvd_isnan(stdl1_ifdf(gp_bool)<0.0006,:) = false;
        logYifc_bprmvd_good_1nan = double(~logYifc_bprmvd_isnan);
        logYifc_bprmvd_good_1nan(logYifc_bprmvd_isnan) = nan;
        %lambda_r_bprmvd_new = 1./mad_bprmvd_rr1.*(Ymodel_bprmvd)./(L_bprmvd*40)*3;
        %lambda_r_bprmvd_new(logYifc_bprmvd_isnan) = 0;
        %lambda_c_bprmvd(logYifc_bprmvd_isnan) = 20;
        %lambda_c_bprmvd([1,L_bprmvd],:) = 0; % no weight for edges
        %logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
        vldpxl = ~any(isnan(logYifc_bprmvd),1);
    otherwise 
        error('not implemented');
end

R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;


Xlogtc_1d = sum(Xlogtc,1);
r_lad = zeros(Ny,L_bprmvd); d_lad = zeros(Ny,L_bprmvd); Rhov_lad = ones(Ny,1);
[logt_est_bprmvd,r_lad(vldpxl,:),d_lad(vldpxl,:),rho_lad,Rhov_lad(vldpxl,:)]...
    = wlad_gadmm_a_v2(Xlogtc_1d(:,vldpxl)', R_bprmvd(:,vldpxl)',...
    'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad);%,...
%     'lambda_r',lambda_r_bprmvd_new(:,vldpxl)');
logt_est_bprmvd = logt_est_bprmvd';
%R_bprmvd_new = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
resNew_bprmvd = R_bprmvd - logt_est_bprmvd*Xlogtc_1d;
% resNewNrm_bprmvd = nansum(nansum(lambda_r_bprmvd_new.*abs(resNew_bprmvd)));
% resNewNrm_bprmvd = nanmedian(lambda_r_bprmvd_new.*abs(resNew_bprmvd).*logYifc_bprmvd_good_1nan,'all');

Ymodel_bprmvd = exp(logAB_bprmvd+logBg_bprmvd+logt_est_bprmvd*Xlogtc_1d);

res_exp = Yifc_bprmvd - Ymodel_bprmvd;
RDimg_bprmvd = if2rd(Ymodel_bprmvd,SFimgc_bprmvd,lbl);
[photon_mad_new] = estimate_photon_noise_CRISM_base(RDimg_bprmvd,wvc_bprmvd,WA_um_pitch_bprmvd,lbl,SFimgc_bprmvd);
mad_bprmvd = robust_v3('med_abs_dev_from_med',res_exp,2,'NOutliers',10);
mad_expected = max(repmat(mad_bprmvd,[1 Ny]),(stdl1_ifdf(gp_bool)+photon_mad_new));
% coeff = 1.75;
% (L_bprmvd/(L_bprmvd-sum(bp_est_bool)))
% lambda_r_bprmvd_new = 1./mad_expected.*(Ymodel_bprmvd)./((L_bprmvd-sum(bp_est_bool))*40);%.*coeff;
lambda_r_bprmvd_new = 1./mad_expected.*(Ymodel_bprmvd)./((L_bprmvd-sum(bp_est_bool))*20);
lambda_r_bprmvd_new(logYifc_bprmvd_isnan) = 0;

resNewNrm_bprmvd = nanmedian(lambda_r_bprmvd_new.*abs(resNew_bprmvd).*logYifc_bprmvd_good_1nan,'all');
%resNrm_bprmvd2 = nanmedian(lambda_r_bprmvd.*abs(resNew_bprmvd).*logYifc_bprmvd_good_1nan_ori,'all');

% lambda_c_bprmvd = 1e-8 * ones(L_bprmvd,Ny);
% lambda_c_bprmvd(logYifc_bprmvd_isnan) = 0.01;
% lambda_c_bprmvd([1,L_bprmvd],:) = 0; % no weight for edges
% edge robust??
% mad_expected_bprmvd_bad = mad_expected./(Ymodel_bprmvd)>0.01;
% lambda_c_bprmvd(2:6,:) = mad_expected_bprmvd_bad(2:6,:)*0.01;
% lambda_c_bprmvd((L_bprmvd-3):(L_bprmvd-1),:) = mad_expected_bprmvd_bad((L_bprmvd-3):(L_bprmvd-1),:)*0.01;

% lambda_c_bprmvd = 1e-8 * ones(L_bprmvd,Ny);
% lambda_c_bprmvd(logYifc_bprmvd_isnan) = 0.01;
% lambda_c_bprmvd([1,L_bprmvd],:) = 0; % no weight for edges
% % edge robust??
% lambda_c_bprmvd([2:4,(L_bprmvd-3):(L_bprmvd-1)],:) = 0.5;


%%
if isdebug
    spcs_bprmvd = logYifc_bprmvd - A_bprmvd(:,idxAlogtc)*Xlogtc;
    spcs_new = logYifc_bprmvd - logt_est_bprmvd*Xlogtc_1d;
    spc_r_bprmvd = logBg_bprmvd + logAB_bprmvd;
    logtc_best = nanmean((A_bprmvd(:,idxAlogtc)*Xlogtc)./sum(Xlogtc,1),2);
    diff_tList(:,1) = logt_est_bprmvd-logtc_best;
    logt_estList = [logtc_best logt_est_bprmvd];
    
    logYifc_bprmvd_bp_1nan = double(isnan(logYifc_bprmvd_good_1nan));
    logYifc_bprmvd_bp_1nan(logYifc_bprmvd_bp_1nan==0) = nan;

    RList(:,1) = vnorms(res_bprmvd,1,2);
    res_nrmList(1) = resNrm_bprmvd;

    figure; ax_tr = subplot(1,1,1); movegui(gcf,'northwest');
    plot(ax_tr,wvc,A(:,idxAlogtc),'Color',[0.5 0.5 0.5],'DisplayName','iter=0'); hold(ax_tr,'on'); 
    figure; ax_dtr = subplot(1,1,1); hold(ax_dtr,'on');movegui(gcf,'north');
    figure; ax_spc = subplot(1,1,1); hold(ax_spc,'off');movegui(gcf,'northeast');
    figure; ax_mdtr = subplot(1,1,1); movegui(gcf,'southwest');
    title(ax_mdtr,'norm of difference of consequtive logt');
    figure; ax_res = subplot(1,1,1); movegui(gcf,'south');
    title(ax_mdtr,'norm of residual');
    figure; ax_resv = subplot(1,1,1); hold(ax_resv,'off');movegui(gcf,'southeast');
    
    bp_est_bool_1nan = double(bp_est_bool);
    bp_est_bool_1nan(bp_est_bool==0) = nan;

    plot(ax_tr,wvc_bprmvd,logt_est_bprmvd,'DisplayName','iter=0');
    for k=422
        %plot(ax_spc,wvc,exp(logYifc_cat(:,k)),'Color','k',...
        %    'DisplayName',sprintf('iter=0;%d cat\n',k));
        %plot(ax_spc,wvc,exp(logYraifc_cat(:,k)),'Color',[0.5 0.5 0.5],...
        %    'DisplayName',sprintf('iter=0;%d cat\n',k));
        hold(ax_spc,'on');
        l1 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)),...
            'DisplayName',sprintf('iter=0;%d\n',k));
        l2 = plot(ax_spc,wvc_bprmvd,exp(spcs_new(:,k)),...
            'DisplayName',sprintf('iter=0;%d\n',k));
        l3 = plot(ax_spc,wvc_bprmvd,exp(spcs_new(:,k)).*logYifc_bprmvd_good_1nan(:,k),'.',...
            'DisplayName',sprintf('iter=0;%d,good\n',k),'Color',l2.Color);
        l3 = plot(ax_spc,wvc_bprmvd,exp(spcs_new(:,k)).*logYifc_bprmvd_bp_1nan(:,k),'x',...
            'DisplayName',sprintf('iter=0;%d,bad\n',k),'Color',l2.Color);
        l4 = plot(ax_spc,wvc_bprmvd,exp(spcs_new(:,k)).*bp_est_bool_1nan,'O',...
            'DisplayName',sprintf('iter=0;%d,bad\n',k),'Color',l2.Color);
        plot(ax_spc,wvc_bprmvd,exp(spc_r_bprmvd(:,k)),':','Color',l1.Color,...
            'DisplayName',sprintf('iter=0;%d m\n',k));
        plot(ax_spc,wvc_bprmvd,exp(logBg_bprmvd(:,k)),'Color',l1.Color,...
            'DisplayName',sprintf('iter=0;%d b\n',k));
        hold off
%             pause

    end
    hold(ax_spc,'off');
    plot(ax_dtr,wvc_bprmvd,diff_tList(:,1),'DisplayName','iter=0');
    plot(ax_res,res_nrmList);
    plot(ax_resv,wvc_bprmvd,res_bprmvd(:,:),'DisplayName','iter=0');
    plot(ax_mdtr,sqrt(mean(diff_tList.^2,1)));

    drawnow;
    fprintf('# bp: % 3d\n',sum(bp_est_bool));
end

%%
% lambda_r_bprmvd = lambda_r_bprmvd_new;
% [ X1(:,vldpxl),Z1(:,vldpxl),C1,~,D1(:,vldpxl),rho,Rhov,~] = huwacbwl1_gadmm_a_v2(...
%     A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
%     'LAMBDA_A',lambda_a_1,'Lambda_c',lambda_c_bprmvd(:,vldpxl),'lambda_r',lambda_r_bprmvd(:,vldpxl),...
%     'tol',1e-5,'maxiter',1000,'verbose',1,'YNormalize',false);
% 
% % lambda_c_bprmvd = 1e-8 * ones(size(C1,2),sum(vldpxl));
%                            
% logBg_bprmvd = C1*Z1;
% % logYifc_model_bprmvd = logBg_bprmvd + A_bprmvd*X1;
% Xlib = X1(idxAlib,:); Xlogtc = X1(idxAlogtc,:);
% logAB_bprmvd = Alib_bprmvd*Xlib;
% R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
% res_bprmvd = R_bprmvd - A_bprmvd(:,idxAlogtc) * Xlogtc;
% %resNrm_bprmvd = nansum(nansum(lambda_r_bprmvd.*abs(res_bprmvd)));
% resNrm_bprmvd = nanmedian(lambda_r_bprmvd.*abs(res_bprmvd).*logYifc_bprmvd_good_1nan,'all');
% 
% 
% switch spike_detection_opt
%     case 'original'
%         Ymodel_bprmvd = exp(logAB_bprmvd+logBg_bprmvd+logtc(gp_bool)*Xlogtc);
%         res_exp1 = Yifc_bprmvd - Ymodel_bprmvd;
%         mad_bprmvd_rr1 = robust_v3('med_abs_dev_from_med',res_exp1,2,'NOutliers',10);
%         mad_bprmvd_log1 = robust_v3('med_abs_dev_from_med',res_bprmvd,2,'NOutliers',10);
%         
%         % First step of de-noising is bad pixel detection.
%         % Our bad pixel detection is based on the median absolute deviation
%         % from the median. If that value is greater than a threshold, then
%         % that will be excluded.
%         % EXCEPTION
%         %  If the residual in the log domain has only a small variance,
%         %  then such pixels will be not considered as a bad pixel. The
%         %  pixel is likely to be the error of a multiplicative components
%         %  in the calibration.
%         logYifc_bprmvd_isnan_bp = false(size(Ymodel_bprmvd));
%         bp_est_bool = and(mad_bprmvd_rr1>0.001,mad_bprmvd_log1>0.005);
%         logYifc_bprmvd_isnan_bp(bp_est_bool,:) = true;
%         logYifc_bprmvd_isnan = or(logYifc_bprmvd_isnan_bp,logYifc_bprmvd_isnan_ori);
%         % for detection of temporal spikes are not evaluated to avoid bias
%         % caused by the initial transmission spectrum or inaccurate BK or
%         % BI.
%         
%         % if the standard deviation of dark frames is small enough, bring back to good one.
% %         logYifc_bprmvd_isnan(stdl1_ifdf(gp_bool)<0.0006,:) = false;
%         logYifc_bprmvd_good_1nan = double(~logYifc_bprmvd_isnan);
%         logYifc_bprmvd_good_1nan(logYifc_bprmvd_isnan) = nan;
%         %lambda_r_bprmvd_new = 1./mad_bprmvd_rr1.*(Ymodel_bprmvd)./(L_bprmvd*40)*3;
%         %lambda_r_bprmvd_new(logYifc_bprmvd_isnan) = 0;
%         lambda_c_bprmvd(logYifc_bprmvd_isnan) = 20;
%         lambda_c_bprmvd([1,L_bprmvd],:) = 0; % no weight for edges
%         %logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
%         vldpxl = ~any(isnan(logYifc_bprmvd),1);
%     otherwise 
%         error('not implemented');
% end
% 
% R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
% 
% 
% Xlogtc_1d = sum(Xlogtc,1);
% r_lad = zeros(Ny,L_bprmvd); d_lad = zeros(Ny,L_bprmvd); Rhov_lad = ones(Ny,1);
% [logt_est_bprmvd,r_lad(vldpxl,:),d_lad(vldpxl,:),rho_lad,Rhov_lad(vldpxl,:)]...
%     = wlad_gadmm_a_v2(Xlogtc_1d(:,vldpxl)', R_bprmvd(:,vldpxl)',...
%     'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad);%,...
% %     'lambda_r',lambda_r_bprmvd_new(:,vldpxl)');
% logt_est_bprmvd = logt_est_bprmvd';
% %R_bprmvd_new = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
% resNew_bprmvd = R_bprmvd - logt_est_bprmvd*Xlogtc_1d;
% % resNewNrm_bprmvd = nansum(nansum(lambda_r_bprmvd_new.*abs(resNew_bprmvd)));
% % resNewNrm_bprmvd = nanmedian(lambda_r_bprmvd_new.*abs(resNew_bprmvd).*logYifc_bprmvd_good_1nan,'all');
% 
% Ymodel_bprmvd = exp(logAB_bprmvd+logBg_bprmvd+logt_est_bprmvd*Xlogtc_1d);
% 
% res_exp = Yifc_bprmvd - Ymodel_bprmvd;
% mad_bprmvd = robust_v3('med_abs_dev_from_med',res_exp,2,'NOutliers',10);
% mad_expected = max(mad_bprmvd,stdl1_ifdf);
% %coeff = nansum(1./stdl1_ifdf)./nansum(1./mad_expected);
% % (L_bprmvd/(L_bprmvd-sum(bp_est_bool)))
% lambda_r_bprmvd_new = 1./mad_expected.*(Ymodel_bprmvd)./((L_bprmvd-sum(bp_est_bool))*40);%.*coeff;
% lambda_r_bprmvd_new(logYifc_bprmvd_isnan) = 0;
% 
% %resNewNrm_bprmvd = nanmedian(lambda_r_bprmvd_new.*abs(resNew_bprmvd).*logYifc_bprmvd_good_1nan,'all');
% %resNrm_bprmvd2 = nanmedian(lambda_r_bprmvd.*abs(resNew_bprmvd).*logYifc_bprmvd_good_1nan_ori,'all');
% if isdebug
%     spcs_bprmvd = logYifc_bprmvd - A_bprmvd(:,idxAlogtc)*Xlogtc;
%     spcs_new = logYifc_bprmvd - logt_est_bprmvd*Xlogtc_1d;
%     spc_r_bprmvd = logBg_bprmvd + logAB_bprmvd;
%     logtc_best = nanmean((A_bprmvd(:,idxAlogtc)*Xlogtc)./sum(Xlogtc,1),2);
%     diff_tList = [diff_tList logt_est_bprmvd-logtc_best];
%     logt_estList = [logtc_best logt_est_bprmvd];
% 
%     RList = [RList vnorms(res_bprmvd,1,2)'];
%     res_nrmList = [res_nrmList resNrm_bprmvd];
% 
%     plot(ax_tr,wvc_bprmvd,logt_est_bprmvd,'DisplayName','iter=0');
%     for k=433
%         %plot(ax_spc,wvc,exp(logYifc_cat(:,k)),'Color','k',...
%         %    'DisplayName',sprintf('iter=0;%d cat\n',k));
%         %plot(ax_spc,wvc,exp(logYraifc_cat(:,k)),'Color',[0.5 0.5 0.5],...
%         %    'DisplayName',sprintf('iter=0;%d cat\n',k));
%         hold(ax_spc,'on');
%         l1 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)),...
%             'DisplayName',sprintf('iter=0;%d\n',k));
%         l2 = plot(ax_spc,wvc_bprmvd,exp(spcs_new(:,k)),...
%             'DisplayName',sprintf('iter=0;%d\n',k));
%         plot(ax_spc,wvc_bprmvd,exp(spc_r_bprmvd(:,k)),':','Color',l1.Color,...
%             'DisplayName',sprintf('iter=0;%d m\n',k));
%         plot(ax_spc,wvc_bprmvd,exp(logBg_bprmvd(:,k)),'Color',l1.Color,...
%             'DisplayName',sprintf('iter=0;%d b\n',k));
%         hold off
% %             pause
% 
%     end
%     hold(ax_spc,'off');
%     plot(ax_dtr,wvc_bprmvd,diff_tList(:,1),'DisplayName','iter=0');
%     plot(ax_res,res_nrmList);
%     plot(ax_resv,wvc_bprmvd,res_bprmvd(:,:),'DisplayName','iter=0');
%     plot(ax_mdtr,sqrt(mean(diff_tList.^2,1)));
% 
%     drawnow;
% end

%% mainloop to improve the estimate
A_bprmvd = [logt_est_bprmvd Alib_bprmvd];
X = [Xlogtc_1d;Xlib];
% D = [sum(D1(idxAlogtc,:),1);D1(idxAlib,:)];
D = [zeros(1,size(D1,2)); D1(idxAlibstrt:end,:)];
C = C1;
Z = Z1;
lambda_a_2 = zeros([size(A_bprmvd,2),1]);

rho = ones([1,Ny]);

% always update lambda_tmp

% lambda_tmp = lambda_tmp*(resNewNrm_bprmvd/resNrm_bprmvd);
%lambda_tmp = lambda_tmp*...
%    nanmedian(lambda_r_bprmvd_new.*logYifc_bprmvd_good_1nan,'all')./nanmedian(lambda_r_bprmvd.*logYifc_bprmvd_good_1nan_ori,'all');

for j=2:nIter+1
    if isdebug
        fprintf('Iter%d,lambda = %3.4e\n',j,lambda_tmp);
    end
    lambda_a_2(2:end) = lambda_tmp;
    lambda_r_bprmvd = lambda_r_bprmvd_new;
    logt_estList(:,j) = A_bprmvd(:,1);
    rr = logYifc_bprmvd - A_bprmvd*X - C*Z;
    if j==2
%         [ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
%             = huwacbl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
%                             'LAMBDA_A',lambda_a_2,'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
%                             'D0',D(:,vldpxl),'X0',X(:,vldpxl),...
%                             'R0',rr(:,vldpxl),...
%                             'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb);
        [ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
            = huwacbwl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
                            'LAMBDA_A',lambda_a_2,'Lambda_R',lambda_r_bprmvd(:,vldpxl),...
                            'Lambda_c',lambda_c_bprmvd(:,vldpxl),'YNormalize',false,...
                            'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
                            'D0',D(:,vldpxl),'X0',X(:,vldpxl),...
                            'R0',rr(:,vldpxl),...
                            'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb);
    else
%        [ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
%            = huwacbl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
%                             'LAMBDA_A',lambda_a_2,'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
%                             'D0',D(:,vldpxl),'X0',X(:,vldpxl),...
%                             'R0',rr(:,vldpxl),'rho',rho(:,vldpxl),'Rhov',Rhov,...
%                             'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb);
       [ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
           = huwacbwl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
                            'LAMBDA_A',lambda_a_2,'Lambda_R',lambda_r_bprmvd(:,vldpxl),...
                            'Lambda_c',lambda_c_bprmvd(:,vldpxl),'YNormalize',false,...
                            'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
                            'D0',D(:,vldpxl),'X0',X(:,vldpxl),...
                            'R0',rr(:,vldpxl),...
                            'rho',rho(:,vldpxl),'Rhov',Rhov,...
                            'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb); 
    end
    
    logBg_bprmvd = C*Z;
    logAB_bprmvd = Alib_bprmvd*X(2:end,:);
    %logYifc_model_bprmvd = logBg_bprmvd + logAB_bprmvd + A_bprmvd(:,1)*X(1,:);
    
    % R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
    res_bprmvd = R_bprmvd-A_bprmvd(:,1)*X(1,:);
    %resNrm_bprmvd = nansum(nansum(lambda_r_bprmvd.*abs(res_bprmvd)));
    resNrm_bprmvd = nanmedian(lambda_r_bprmvd.*abs(res_bprmvd).*logYifc_bprmvd_good_1nan,'all');
    
    RDimg_bprmvd = if2rd(Ymodel_bprmvd,SFimgc_bprmvd,lbl);
    [photon_mad_new] = estimate_photon_noise_CRISM_base(RDimg_bprmvd,wvc_bprmvd,WA_um_pitch_bprmvd,lbl,SFimgc_bprmvd);

    
    switch spike_detection_opt
        case 'original'
            Ymodel_bprmvd = exp(logAB_bprmvd+logBg_bprmvd+A_bprmvd(:,1)*X(1,:));
            res_exp = Yifc_bprmvd - Ymodel_bprmvd;
    %         rr1 = logYifc_bprmvd_ori-logYifc_model_bprmvd;
            mad_bprmvd_rr = robust_v3('med_abs_dev_from_med',res_exp,2,'NOutliers',10);
            mad_expected = max(repmat(mad_bprmvd,[1 Ny]),(stdl1_ifdf(gp_bool)+photon_mad_new));
            % mad_expected = max(mad_bprmvd_rr,stdl1_ifdf);
            % rr1_std = nanstd(rr1,[],2);
            %logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(res_exp1)>0.01 );

            % More rigid bad pixel detection
            logYifc_bprmvd_isnan_bp = false(size(Ymodel_bprmvd));
            bp_est_bool = mad_bprmvd_rr>0.001;
            logYifc_bprmvd_isnan_bp(bp_est_bool,:) = true;
            
            % Now perform temporal spike removal (assuming that bias problem is gone)
            res_exp_scaled = res_exp./(exp(A_bprmvd(:,1)).^X(1,:));
            logYifc_bprmvd_isnan_spk = abs(res_exp)>0.0015;
            
            logYifc_bprmvd_isnan = or(or(logYifc_bprmvd_isnan_bp,logYifc_bprmvd_isnan_ori),...
                logYifc_bprmvd_isnan_spk);

            % if the standard deviation of dark frames is small enough, bring back to good one.
    %         logYifc_bprmvd_isnan(stdl1_ifdf(gp_bool)<0.0006,:) = false;
            logYifc_bprmvd_good_1nan = double(~logYifc_bprmvd_isnan);
            logYifc_bprmvd_good_1nan(logYifc_bprmvd_isnan) = nan;
            lambda_r_bprmvd_new = 1./mad_expected.*(Ymodel_bprmvd)./((L_bprmvd-sum(bp_est_bool))*20);
            lambda_r_bprmvd_new(logYifc_bprmvd_isnan) = 0;
            
            %lambda_c_bprmvd = 1e-8 * ones(L_bprmvd,Ny);
            %lambda_c_bprmvd(logYifc_bprmvd_isnan) = 0.01;
            %lambda_c_bprmvd([1,L_bprmvd],:) = 0; % no weight for edges
            % edge robust??
            %lambda_c_bprmvd([2:4,(L_bprmvd-3):(L_bprmvd-1)],:) = 0.5;
            
             %lambda_c_bprmvd = 1e-8 * ones(L_bprmvd,Ny);
             %lambda_c_bprmvd(logYifc_bprmvd_isnan) = 0.01;
             %lambda_c_bprmvd([1,L_bprmvd],:) = 0; % no weight for edges
            % edge robust??
%             mad_expected_bprmvd_bad = mad_expected./(Ymodel_bprmvd)>0.01;
%             lambda_c_bprmvd(2:6,:) = mad_expected_bprmvd_bad(2:6,:)*0.01;
%             lambda_c_bprmvd((L_bprmvd-3):(L_bprmvd-1),:) = mad_expected_bprmvd_bad((L_bprmvd-3):(L_bprmvd-1),:)*0.01;
            
            % lambda_r_bprmvd_new(logYifc_bprmvd_isnan_ori) = 0;
            % lambda_r_bprmvd_new(logYifc_bprmvd_isnan_spk) = 0;
            % lambda_c_bprmvd(logYifc_bprmvd_isnan) = 20;
            % lambda_c_bprmvd([1,L_bprmvd],:) = 0; % no weight for edges
            %logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
            vldpxl = (sum(logYifc_bprmvd_isnan,1)/L_bprmvd) < 0.7;
        otherwise 
            error('not implemented');
    end
    
    
%     switch spike_detection_opt
%         case 'original'
%             rr = logYifc_bprmvd-logYifc_model_bprmvd;
%             logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(rr)>0.015 );
%             lambda_r_bprmvd_new = lambda_r_bprmvd;
%             lambda_r_bprmvd_new(logYifc_bprmvd_isnan) = 0;
%             lambda_c_bprmvd(logYifc_bprmvd_isnan) = 20;
%             lambda_c_bprmvd([1,L_bprmvd],:) = 0; % no weight for edges
%             % logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
%             vldpxl = (sum(logYifc_bprmvd_isnan,1)/L_bprmvd) < 0.8;
%         otherwise 
%             error('not implemented');
%     end   

    % remove nans
    
    
    %update logt_est!
    R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
%     [logt_est_bprmvd,r_lad(vldpxl,:),d_lad(vldpxl,:),rho_lad,Rhov_lad(vldpxl,:)]...
%         = wlad_gadmm_a_v2(X(1,vldpxl)', R_bprmvd(:,vldpxl)',...
%            'lambda_r',lambda_r_bprmvd_new(:,vldpxl)',...
%            'X0',A_bprmvd(:,1)','D0',d_lad(vldpxl,:),...
%            'rho',rho_lad,'Rhov',Rhov_lad(vldpxl,:),...
%            'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad);
    
    [logt_est_bprmvd,r_lad(vldpxl,:),d_lad(vldpxl,:),rho_lad,Rhov_lad(vldpxl,:)]...
    = wlad_gadmm_a_v2(X(1,vldpxl)', R_bprmvd(:,vldpxl)','lambda_r',lambda_r_bprmvd_new(:,vldpxl)',...
    'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad);
       
    logt_est_bprmvd = logt_est_bprmvd';
%     logt_est = update_logt_est(R,X(1,:));
    diff_tList(:,j) = logt_est_bprmvd-logt_estList(:,j);
    
    if (sqrt(mean(diff_tList(:,j).^2,1))/sqrt(mean(diff_tList(:,1).^2,1)) < 0.05)
        % break;
    end
    
    
    %R_bprmvdNew = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
    resNew_bprmvd = R_bprmvd - logt_est_bprmvd*X(1,:);
    %resNewNrm_bprmvd = nansum(nansum(lambda_r_bprmvd_new.*abs(resNew_bprmvd)));
    resNewNrm_bprmvd = nanmedian(lambda_r_bprmvd_new.*abs(resNew_bprmvd).*logYifc_bprmvd_good_1nan,'all');
    
    A_bprmvd(:,1) = logt_est_bprmvd;
    
    %lambda_tmp = lambda_tmp*resNewNrm_bprmvd/resNrm_bprmvd;

        
    if isdebug
        res_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd-A_bprmvd(:,1)*X(1,:);
        % resNrm_bprmvd = nansum(nansum(res_bprmvd.^2));
        RList =[RList vnorms(res_bprmvd,1,2)'];
        res_nrmList= [res_nrmList resNewNrm_bprmvd];
        logt_estList = [logt_estList logt_est_bprmvd];
        diff_tList = [diff_tList logt_est_bprmvd-logt_estList(:,end-1)];
        spcs_bprmvd = logYifc_bprmvd - A_bprmvd(:,1)*X(1,:);
        spc_r_bprmvd = logBg_bprmvd + logAB_bprmvd;
        
        logYifc_bprmvd_bp_1nan = double(isnan(logYifc_bprmvd_good_1nan));
        logYifc_bprmvd_bp_1nan(logYifc_bprmvd_bp_1nan==0) = nan;
        
        bp_est_bool_1nan = double(bp_est_bool);
        bp_est_bool_1nan(bp_est_bool==0) = nan;
        
        plot(ax_tr,wvc_bprmvd,logt_est_bprmvd,'DisplayName',sprintf('iter=%d',j));
        for k=422
            hold(ax_spc,'on');
            %plot(ax_spc,wvc,exp(logYifc_cat(:,k)),'Color','k',...
            %    'DisplayName',sprintf('iter=0;%d cat\n',k));
            hold(ax_spc,'on');
            l1 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)),...
                            'DisplayName',sprintf('iter=%d;%d\n',j,k));
            l3 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)).*logYifc_bprmvd_good_1nan(:,k),'.',...
                'DisplayName',sprintf('iter=0;%d\n',k),'Color',l1.Color);
            l4 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)).*logYifc_bprmvd_bp_1nan(:,k),'x',...
                'DisplayName',sprintf('iter=0;%d,good\n',k),'Color',l1.Color);
            plot(ax_spc,wvc_bprmvd,exp(spc_r_bprmvd(:,k)),':','Color',l1.Color,...
                            'DisplayName',sprintf('iter=0;%d m\n',j,k));
            plot(ax_spc,wvc_bprmvd,exp(logBg_bprmvd(:,k)),'--','Color',l1.Color,...
                'DisplayName',sprintf('iter=%d;%d b\n',j,k));
            l4 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)).*bp_est_bool_1nan,'O',...
            'DisplayName',sprintf('iter=0;%d,bad\n',k),'Color',l1.Color);
        end
        hold(ax_spc,'off');
        plot(ax_dtr,wvc_bprmvd,diff_tList(:,j),'DisplayName',sprintf('iter=%d',j));
        plot(ax_res,res_nrmList);
        plot(ax_resv,wvc_bprmvd,res_bprmvd(:,:),'DisplayName',sprintf('iter=%d',j));
        plot(ax_mdtr,sqrt(mean(diff_tList.^2,1)));
        drawnow;
        fprintf('# bp: % 3d\n',sum(bp_est_bool));
    end
end

res_exp_scaled = res_exp./(exp(A_bprmvd(:,1)).^X(1,:));
logYifc_bprmvd_isnan_spk = abs(res_exp_scaled)>0.0015;
logYifc_bprmvd_isnan = or(or(logYifc_bprmvd_isnan_bp,logYifc_bprmvd_isnan_ori),...
                logYifc_bprmvd_isnan_spk);
logYifc_bprmvd_good_1nan = double(~logYifc_bprmvd_isnan);
logYifc_bprmvd_good_1nan(logYifc_bprmvd_isnan) = nan;

% lambda_c_bprmvd(logYifc_bprmvd_isnan) = 0.01;
%% last iteration after estimating log_est
rr = logYifc_bprmvd - A_bprmvd*X - C*Z;
lambda_a_2(2:end) = lambda_tmp;
lambda_r_bprmvd = lambda_r_bprmvd_new;

[ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
   = huwacbwl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
                    'LAMBDA_A',lambda_a_2,'Lambda_R',lambda_r_bprmvd(:,vldpxl),...
                    'Lambda_c',lambda_c_bprmvd(:,vldpxl),'YNormalize',false,...
                    'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
                    'D0',D(:,vldpxl),'X0',X(:,vldpxl),...
                    'R0',rr(:,vldpxl),...
                    'rho',rho(:,vldpxl),'Rhov',Rhov,...
                    'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb); 

% [ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
%     = huwacbl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
%                             'LAMBDA_A',lambda_a_2,'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
%                             'D0',D(:,vldpxl),'X0',X(:,vldpxl),'R0',rr(:,vldpxl),...
%                             'rho',rho(:,vldpxl),'Rhov',Rhov,...
%                             'verbose',verbose_huwacb,...
%                             'tol',tol_huwacb,'maxiter',maxiter_huwacb);

% original wavelength channels
logBg = nan(size(logYifc));
logBg(gp_bool,:) = C*Z;

logAB = Alib*X(2:end,:);

logt_est = nan(size(logYifc,1),1);
logt_est(gp_bool) = A_bprmvd(:,1);

if isdebug
    logBg_bprmvd = C*Z;
    logAB_bprmvd = Alib_bprmvd*X(2:end,:);
    res_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd-A_bprmvd(:,1)*X(1,:);
    % resNrm_bprmvd = nansum(nansum(res_bprmvd.^2));
    RList =[RList vnorms(res_bprmvd,1,2)'];
    res_nrmList= [res_nrmList resNewNrm_bprmvd];
%     logt_estList = [logt_estList logt_est_bprmvd];
%     diff_tList = [diff_tList logt_est_bprmvd-logt_estList(:,end-1)];
    spcs_bprmvd = logYifc_bprmvd - A_bprmvd(:,1)*X(1,:);
    spc_r_bprmvd = logBg_bprmvd + logAB_bprmvd;

    logYifc_bprmvd_bp_1nan = double(isnan(logYifc_bprmvd_good_1nan));
    logYifc_bprmvd_bp_1nan(logYifc_bprmvd_bp_1nan==0) = nan;
    
    bp_est_bool_1nan = double(bp_est_bool);
    bp_est_bool_1nan(bp_est_bool==0) = nan;

%     plot(ax_tr,wvc_bprmvd,logt_est_bprmvd,'DisplayName',sprintf('iter=%d',j));
    for k=422
        hold(ax_spc,'on');
        %plot(ax_spc,wvc,exp(logYifc_cat(:,k)),'Color','k',...
        %    'DisplayName',sprintf('iter=0;%d cat\n',k));
        hold(ax_spc,'on');
        l1 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)),...
                        'DisplayName',sprintf('iter=%d;%d\n',j,k));
        l3 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)).*logYifc_bprmvd_good_1nan(:,k),'.',...
            'DisplayName',sprintf('iter=0;%d\n',k),'Color',l1.Color);
        l4 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)).*logYifc_bprmvd_bp_1nan(:,k),'x',...
            'DisplayName',sprintf('iter=0;%d,good\n',k),'Color',l1.Color);
        plot(ax_spc,wvc_bprmvd,exp(spc_r_bprmvd(:,k)),':','Color',l1.Color,...
                        'DisplayName',sprintf('iter=0;%d m\n',j,k));
        plot(ax_spc,wvc_bprmvd,exp(logBg_bprmvd(:,k)),'--','Color',l1.Color,...
            'DisplayName',sprintf('iter=%d;%d b\n',j,k));
        l4 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)).*bp_est_bool_1nan,'O',...
            'DisplayName',sprintf('iter=0;%d,bad\n',k),'Color',l1.Color);
    end
    hold(ax_spc,'off');
%     plot(ax_dtr,wvc_bprmvd,diff_tList(:,j),'DisplayName',sprintf('iter=%d',j));
%     plot(ax_res,res_nrmList);
%     plot(ax_resv,wvc_bprmvd,res_bprmvd(:,:),'DisplayName',sprintf('iter=%d',j));
%     plot(ax_mdtr,sqrt(mean(diff_tList.^2,1)));
    drawnow;
end



% corrected spectra
logYifc_cor = nan(size(logYifc));
logYifc_bprmvd(logYifc_bprmvd_isnan) = nan;
logYifc_cor(gp_bool,:) = logYifc_bprmvd - A_bprmvd(:,1)*X(1,:);
% corrected spectra
logYifc_cor_ori = nan(size(logYifc));
logYifc_cor_ori(gp_bool,:) = logYifc_bprmvd_ori - A_bprmvd(:,1)*X(1,:);

logYifc_isnan = true(size(logYifc));
logYifc_isnan(gp_bool,:) = logYifc_bprmvd_isnan;

bp_bool = ~gp_bool;
logBg = interp_nan_column(logBg,repmat(bp_bool,[1,Ny]),wvc);

% residual
rr_ori = nan(size(logYifc));
rr_ori(gp_bool,:) = logYifc_bprmvd_ori - logAB(gp_bool,:) - logBg(gp_bool,:) - A_bprmvd(:,1)*X(1,:);

vldpxl = (sum(logYifc_bprmvd_isnan,1)/L_bprmvd) < 0.7;

Ymodel = exp(logAB+logBg+logt_est*X(1,:));

res_exp = Yif - Ymodel;

std_expected_scaled = nan(size(logYifc));
std_expected_scaled(gp_bool,:) = mad_expected ./ exp(A_bprmvd(:,1).*X(1,:)) ./ norminv(0.75);

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

