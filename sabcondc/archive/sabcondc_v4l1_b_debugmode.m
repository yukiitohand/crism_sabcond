% Aice will be kept until the end
function [ logt_est,logYifc_cor,logAB,logBg,logIce,logYifc_cor_ori,logYifc_isnan,ancillary,rr_ori,vldpxl] = sabcondc_v4l1_b_debugmode( Alib,Aice,logYifc,wvc,logtc,varargin )
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


if isdebug
    wvc_plot = nan(size(wvc));
    wvc_plot(gp_bool) = wvc_bprmvd;
    logtc_plot = nan(size(A(:,idxAlogtc)));
    logtc_plot(gp_bool,:) = A_bprmvd(:,idxAlogtc);
    fig_t = figure;
    plot(wvc_plot,logtc_plot,'.-');
    set(gca,'FontSize',18);
    xlim([1000 2650]);
    xlabel('Wavelength [nm]');
    ax_t = gca;
    drawnow;
    export_fig(fig_t,'Iter0_t_est.png','-r150','-transparent');
end

% logBg_bprmvd = nan(size(logYifc));
%% initialization of logt_est
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

Xlib = X1(idxAlib,:); Xlogtc = X1(idxAlogtc,:);
Xice = X1(idxAice,:);
logAB_bprmvd = Alib_bprmvd*Xlib;
logIce_bprmvd = Aice_bprmvd*Xice;

if isdebug
    spc_cor = logYifc_bprmvd - logIce_bprmvd - A_bprmvd(:,idxAlogtc) * Xlogtc;
    spc_cor_plot = nan(size(logYifc));
    spc_cor_plot(gp_bool,:) = spc_cor;
    fig_spc = figure;
    plot(wvc_plot,spc_cor_plot,'.-');
    xlim([1000 2650]);
    set(gca,'FontSize',18);
    ax_spc = gca;
    xlabel('Wavelength [nm]');
    drawnow;
    export_fig(fig_spc,'Iter1_spc_cor_iter.png','-r150','-transparent');
    
    spc_model = logAB_bprmvd+logBg_bprmvd;
    spc_model_plot = nan(size(logYifc));
    spc_model_plot(gp_bool,:) = spc_model;
    fig_spc_model = figure;
    plot(wvc_plot,spc_model_plot,'.-');
    xlim([1000 2650]);
    set(gca,'FontSize',18);
    ax_spc_model = gca;
    xlabel('Wavelength [nm]');
    drawnow;
    export_fig(fig_spc_model,'Iter1_spc_model_iter.png','-r150','-transparent');
    
    spc_Ice_plot = nan(size(logYifc));
    spc_Ice_plot(gp_bool,:) = logIce_bprmvd;
    fig_ice = figure;
    plot(wvc_plot,spc_Ice_plot,'.-');
    xlim([1000 2650]);
    set(gca,'FontSize',18);
    ax_Ice_model = gca;
    set(ax_Ice_model,'YLim',ax_t.YLim);
    xlabel('Wavelength [nm]');
    drawnow;
    export_fig(fig_ice,'Iter1_spc_Ice.png','-r150','-transparent');
    
end


switch spike_detection_opt
%     case 'incremental'
%         rr1 = logYifc_bprmvd-logYifc_model_bprmvd;
%         logYifc_bprmvd_isnan = or(logYifc_bprmvd_isnan,abs(rr1)>0.05);
%         logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
    case 'original'
        rr1 = logYifc_bprmvd_ori-logYifc_model_bprmvd;
        if isdebug
            fig_residual = figure;
            rr1_plot = nan(size(logYifc));
            rr1_plot(gp_bool,:) = rr1;
            plot(wvc_plot,rr1_plot,'.-','Color',[0.5 0.5 0.5]);
            xlim([1000 2650]);
            set(gca,'FontSize',18);
            ax_residual = gca;
            xlabel('Wavelength [nm]');
            drawnow;
            export_fig(fig_residual,'Iter1_residual.png','-r150','-transparent');

            hold(ax_residual,'on');
        end
        rr1_std = std(rr1,[],2);
        logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(rr1)>0.1 );
        logYifc_bprmvd_isnan(rr1_std<0.015,:) = false; % if the standard deviation of the channel is small enough, bring back to good one.
        logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
        vldpxl = ~any(isnan(logYifc_bprmvd),1);
        if isdebug
            logYifc_bprmvd_isnan_tmp = double(logYifc_bprmvd_isnan==1);
            logYifc_bprmvd_isnan_tmp(logYifc_bprmvd_isnan_tmp==0) = nan;
            logYifc_bprmvd_isnan_tmp_plot = nan(size(logYifc));
            logYifc_bprmvd_isnan_tmp_plot(gp_bool,:) = logYifc_bprmvd_isnan_tmp;
            plot(ax_residual,wvc_plot, rr1_plot.*logYifc_bprmvd_isnan_tmp_plot,'X','Color','r');
            xlim(ax_residual,[1000 2650]);
            set(ax_residual,'FontSize',18);
            xlabel(ax_residual,'Wavelength [nm]');
            drawnow;
            export_fig(fig_residual,'Iter1_residual_detect.png','-r150','-transparent');
            hold(ax_residual,'off');
        end
%     case 'so_original'
%         rr1 = logYifc - logYifc_model_bprmvd;
%         logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(rr1)>0.05 );
    otherwise 
        error('not implemented');
end


R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd - logIce_bprmvd;
res_bprmvd = R_bprmvd - A_bprmvd(:,idxAlogtc) * Xlogtc;
resNrm_bprmvd = nansum(nansum(abs(res_bprmvd)));

if isdebug
    spc_cor = logYifc_bprmvd - logIce_bprmvd - A_bprmvd(:,idxAlogtc) * Xlogtc;
    spc_cor_plot(gp_bool,:) = spc_cor;
    plot(ax_spc,wvc_plot,spc_cor_plot,'.-');
    xlim(ax_spc,[1000 2650]);
    set(ax_spc,'FontSize',18);
    xlabel(ax_spc,'Wavelength [nm]');
    drawnow;
    export_fig(fig_spc,'Iter1_spc_cor_denoised.png','-r150','-transparent');
end


Xlogtc_1d = sum(Xlogtc,1);
r_lad = nan(Ny,L_bprmvd); d_lad = nan(Ny,L_bprmvd); Rhov_lad = ones(Ny,1);
[logt_est_bprmvd,r_lad(vldpxl,:),d_lad(vldpxl,:),rho_lad,Rhov_lad(vldpxl,:)] = lad_gadmm_a_v2(Xlogtc_1d(:,vldpxl)', R_bprmvd(:,vldpxl)',...
                                'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad);
logt_est_bprmvd = logt_est_bprmvd';

resNew_bprmvd = R_bprmvd - logt_est_bprmvd*Xlogtc_1d;
resNewNrm_bprmvd = nansum(nansum(abs(resNew_bprmvd)));

if isdebug
    logt_est_plot = nan([size(logYifc,1),1]);
    logt_est_plot(gp_bool) = logt_est_bprmvd;
    plot(ax_t,wvc_plot,logt_est_plot,'r.-');
    xlim(ax_t,[1000 2650]);
    set(ax_t,'FontSize',18);
    xlabel(ax_t,'Wavelength [nm]');
    drawnow;
    export_fig(fig_t,'Iter1_t_est.png','-r150','-transparent');
end


%% mainloop to improve the estimate
A_bprmvd = [logt_est_bprmvd Aice_bprmvd Alib_bprmvd];
nA = size(A_bprmvd,2);
X = [Xlogtc_1d;Xice;Xlib];
D = [zeros(1,size(D1,2));D1(idxAicestrt:end,:)];

C = C1;
Z = Z1;
lambda_a_2 = zeros([1,nA]);

idxAice2 = false(1,nA); idxAice2(2:(nIce+1)) = true;
idxAlib2 = false(1,nA); idxAlib2((nIce+2):end) = true;

rho = ones([1,Ny]);

% always update lambda_tmp
lambda_tmp = lambda_tmp*resNewNrm_bprmvd/resNrm_bprmvd;
lam_tmp_ice = lam_tmp_ice*resNewNrm_bprmvd/resNrm_bprmvd;

for j=2:nIter
    if isdebug
        fprintf('Iter%d,lambda = %3.4e\n',j,lambda_tmp);
    end
    
    if ~isempty(Aice)
        lambda_a_2(idxAlib2) = lambda_tmp;
        lambda_a_2(idxAice2) = lam_tmp_ice;
    else
        lambda_a_2(idxAlib2) = lambda_tmp;  
    end 

    logt_estList(:,j) = A_bprmvd(:,1);
    rr = logYifc_bprmvd - A_bprmvd*X - C*Z;
    if j==2
        [ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
            = huwacbl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
                            'LAMBDA_A',lambda_a_2,'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
                            'D0',D(:,vldpxl),'X0',X(:,vldpxl),...
                            'R0',rr(:,vldpxl),...
                            'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb);
    else
       [ X(:,vldpxl),Z(:,vldpxl),C,~,D(:,vldpxl),rho(:,vldpxl),Rhov ]...
           = huwacbl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,...
                            'LAMBDA_A',lambda_a_2,'CONCAVEBASE',C,'Z0',Z(:,vldpxl),...
                            'D0',D(:,vldpxl),'X0',X(:,vldpxl),...
                            'R0',rr(:,vldpxl),'rho',rho(:,vldpxl),'Rhov',Rhov,...
                            'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb); 
    end
    
    logBg_bprmvd = C*Z;
    logAB_bprmvd = Alib_bprmvd*X(idxAlib2,:);
    logIce_bprmvd = Aice_bprmvd*X(idxAice2,:);
    logYifc_model_bprmvd = logBg_bprmvd + logAB_bprmvd + A_bprmvd(:,1)*X(1,:) + logIce_bprmvd;
    
    if isdebug
        spc_cor = logYifc_bprmvd - logIce_bprmvd - A_bprmvd(:,1) * X(1,:);
        spc_cor_plot(gp_bool,:) = spc_cor;
        plot(ax_spc,wvc_plot,spc_cor_plot,'.-');
        xlim(ax_spc,[1000 2650]);
        set(ax_spc,'FontSize',18);
        xlabel(ax_spc,'Wavelength [nm]');
        drawnow;
        export_fig(fig_spc,sprintf('Iter%d_spc_cor.png',j),'-r150','-transparent');
        
        spc_model = logAB_bprmvd+logBg_bprmvd;
        spc_model_plot(gp_bool,:) = spc_model;
        plot(ax_spc_model,wvc_plot,spc_model_plot,'.-');
        xlim(ax_spc_model,[1000 2650]);
        set(ax_spc_model,'FontSize',18);
        xlabel(ax_spc_model,'Wavelength [nm]');
        drawnow;
        export_fig(fig_spc_model,sprintf('Iter%d_spc_model.png',j),'-r150','-transparent');
        
        spc_Ice_plot = nan(size(logYifc));
        spc_Ice_plot(gp_bool,:) = logIce_bprmvd;
        plot(ax_Ice_model,wvc_plot,spc_Ice_plot,'.-');
        xlim(ax_Ice_model,[1000 2650]);
        set(ax_Ice_model,'FontSize',18);
        set(ax_Ice_model,'YLim',ax_t.YLim);
        xlabel(ax_Ice_model,'Wavelength [nm]');
        drawnow;
        export_fig(fig_ice,sprintf('Iter%d_spc_Ice.png',j),'-r150','-transparent');
    

    end
    
    switch spike_detection_opt
% %         case 'incremental'
% %             rr = logYifc_bprmvd-logYifc_model_bprmvd;
% %             logYifc_bprmvd_isnan = or(logYifc_bprmvd_isnan,abs(rr)>0.015);
% %             logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
        case 'original'
            rr = logYifc_bprmvd_ori-logYifc_model_bprmvd;
            if isdebug
                rr_tmp = logYifc_bprmvd-logYifc_model_bprmvd;
                rr_tmp_plot = nan(size(logYifc));
                rr_tmp_plot(gp_bool,:) = rr_tmp;
                plot(ax_residual,wvc_plot,rr_tmp_plot,'.-','Color',[0.5 0.5 0.5]);
                xlim(ax_residual,[1000 2650]);
                set(ax_residual,'FontSize',18);
                hold(ax_residual,'on');
                xlabel(ax_residual,'Wavelength [nm]');
                drawnow;
                export_fig(fig_residual,sprintf('Iter%d_residual.png',j),'-r150','-transparent');

                logYifc_bprmvd_isnan_old = logYifc_bprmvd_isnan;
            end
            logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(rr)>0.015 );
            logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
            vldpxl = (sum(logYifc_bprmvd_isnan,1)/L_bprmvd) < 0.8;
            if isdebug
                logYifc_bprmvd_isnan_tmp = double(and(logYifc_bprmvd_isnan==1,logYifc_bprmvd_isnan_old==0));
                logYifc_bprmvd_isnan_tmp(logYifc_bprmvd_isnan_tmp==0) = nan;
                logYifc_bprmvd_isnan_tmp_plot(gp_bool,:) = logYifc_bprmvd_isnan_tmp;
                plot(ax_residual,wvc_plot, rr_tmp_plot.*logYifc_bprmvd_isnan_tmp_plot,'X','Color','r');
                xlim(ax_residual,[1000 2650]);
                set(ax_residual,'FontSize',18);
                xlabel(ax_residual,'Wavelength [nm]');
                drawnow;
                export_fig(fig_residual,sprintf('Iter%d_residual_detect.png',j),'-r150','-transparent');
                hold(ax_residual,'off');
            end
        otherwise 
            error('not implemented');
    end
    
    if isdebug
        spc_cor = logYifc_bprmvd - logIce_bprmvd - A_bprmvd(:,1) * X(1,:);
        spc_cor_plot(gp_bool,:) = spc_cor;
        plot(ax_spc,wvc_plot,spc_cor_plot,'.-');
        xlim(ax_spc,[1000 2650]);
        set(ax_spc,'FontSize',18);
        xlabel(ax_spc,'Wavelength [nm]');
        drawnow;
        export_fig(fig_spc,sprintf('Iter%d_spc_cor_denoised.png',j),'-r150','-transparent');
    end
    
    %update logt_est!
    R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd - logIce_bprmvd;
    [logt_est_bprmvd,r_lad(vldpxl,:),d_lad(vldpxl,:),rho_lad,Rhov_lad(vldpxl,:)]...
        = lad_gadmm_a_v2(X(1,vldpxl)', R_bprmvd(:,vldpxl)',...
        'X0',A_bprmvd(:,1)','R0',r_lad(vldpxl,:),'D0',d_lad(vldpxl,:),...
        'rho',rho_lad,'Rhov',Rhov_lad(vldpxl,:),...
        'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad);
    logt_est_bprmvd = logt_est_bprmvd';
    
%     logt_est_ice = lad_gadmm_a_v2(X(1,ii(1:100))', logIce_bprmvd(:,vldpxl(ii(1:100)))',...
%                                'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad);
%     logt_est_ice = logt_est_ice';
%     logt_est_ice = logIce_bprmvd(:,vldpxl(ii(1:30))) * X(1,ii(1:30))' / (X(1,ii(1:30))*X(1,ii(1:30))');
%     logt_est_bprmvd = logt_est_bprmvd + logt_est_ice;
    
%     logt_est = update_logt_est(R,X(1,:));
    diff_tList(:,j) = logt_est_bprmvd-logt_estList(:,j);
    
    if (sqrt(mean(diff_tList(:,j).^2,1))/sqrt(mean(diff_tList(:,1).^2,1)) < 0.05)
        % break;
    end
    
    res_bprmvd = R_bprmvd-A_bprmvd(:,1)*X(1,:);
    resNrm_bprmvd = nansum(nansum(abs(res_bprmvd)));
    resNew_bprmvd = R_bprmvd - logt_est_bprmvd*X(1,:);
    resNewNrm_bprmvd = nansum(nansum(abs(resNew_bprmvd)));
    
    A_bprmvd(:,1) = logt_est_bprmvd;
    
    lambda_tmp = lambda_tmp*resNewNrm_bprmvd/resNrm_bprmvd;
    lam_tmp_ice = lam_tmp_ice*resNewNrm_bprmvd/resNrm_bprmvd;

    if isdebug
        logt_est_plot(gp_bool) = logt_est_bprmvd;
        plot(ax_t,wvc_plot,logt_est_plot,'r.-');
        xlim(ax_t,[1000 2650]);
        set(ax_t,'FontSize',18);
        xlabel(ax_t,'Wavelength [nm]');
        drawnow;
        export_fig(fig_t,sprintf('Iter%d_t_est.png',j),'-r150','-transparent');
    end


end

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

logAB = Alib*X(idxAlib2,:);
logIce = Aice*X(idxAice2,:);

logt_est = nan(size(logYifc,1),1);
logt_est(gp_bool) = A_bprmvd(:,1);

% corrected spectra
logYifc_cor = nan(size(logYifc));
logYifc_bprmvd(logYifc_bprmvd_isnan) = nan;
logYifc_cor(gp_bool,:) = logYifc_bprmvd - A_bprmvd(:,1)*X(1,:) - logIce(gp_bool,:);
% corrected spectra
logYifc_cor_ori = nan(size(logYifc));
logYifc_cor_ori(gp_bool,:) = logYifc_bprmvd_ori - A_bprmvd(:,1)*X(1,:) - logIce(gp_bool,:);

logYifc_isnan = true(size(logYifc));
logYifc_isnan(gp_bool,:) = logYifc_bprmvd_isnan;

logBg = interp_nan_column(logBg,logYifc_isnan,wvc);

% residual
rr_ori = nan(size(logYifc));
rr_ori(gp_bool,:) = logYifc_bprmvd_ori - logAB(gp_bool,:) - logBg(gp_bool,:) - logIce(gp_bool,:) - A_bprmvd(:,1)*X(1,:);

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