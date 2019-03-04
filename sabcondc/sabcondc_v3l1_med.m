function [ logt_est,logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,ancillary,rr_ori,vldpxl] = sabcondc_v3l1_med( Alib,logYifc,wvc,logtc,varargin )
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

A = [logtc Alib];
idxAlibstrt = size(logtc,2)+1;
idxAlogtc = false(1,size(A,2));
idxAlogtc(1:idxAlibstrt-1) = true;
idxAlib = false(1,size(A,2));
idxAlib(idxAlibstrt:end) = true;

N_A1 = size(A,2);


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

% logBg_bprmvd = nan(size(logYifc));
%% initialization of logt_est
lambda_a_1 = zeros([1,size(A,2)]);
lambda_tmp = lambda_a;
% lambda_a_1(idxAlibstrt:end) = lambda_tmp;
lambda_a_1(idxAlib) = lambda_tmp;

X1 = zeros(N_A1,Ny);
Z1 = zeros(L_bprmvd,Ny);
D1 = zeros([N_A1+L_bprmvd*2,Ny]);

[ X1(:,vldpxl),Z1(:,vldpxl),C1,~,D1(:,vldpxl),rho,Rhov,~] = huwacbl1_gadmm_a_v2(A_bprmvd,logYifc_bprmvd(:,vldpxl),wvc_bprmvd,'LAMBDA_A',lambda_a_1,...
                               'tol',1e-5,'maxiter',maxiter_huwacb,'verbose',verbose_huwacb);
                           
logBg_bprmvd = C1*Z1;
logYifc_model_bprmvd = logBg_bprmvd + A_bprmvd*X1;
Xlib = X1(idxAlib,:); Xlogtc = X1(idxAlogtc,:);
logAB_bprmvd = Alib_bprmvd*Xlib;
R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
res_bprmvd = R_bprmvd - A_bprmvd(:,idxAlogtc) * Xlogtc;
%resNrm_bprmvd = nansum(nansum(abs(res_bprmvd)));
resNrm_bprmvd = nanmedian(nanmedian(abs(res_bprmvd)));

% lambdaList = [1e-3 2e-3 4e-3 8e-3 1e-2];
% lambdaList = [1e-5 1e-4 1e-3 1e-2 1e-1];
% figure(1); plot(logYifc_bprmvd(:,102));
% for i=1:length(lambdaList)
%     lambda_a_1(idxAlibstrt:end) = lambdaList(i);
%     [ X1,Z1,C1,rtmp,~ ] = huwacbl1_gadmm_a(A_bprmvd,logYifc_bprmvd,wvc_bprmvd,'LAMBDA_A',lambda_a_1);
%     logBg = C1*Z1;
%     figure(1);
%     hold on;
%     plot(logBg(:,102),'DisplayName',num2str(lambdaList(i)));
%     figure(2);
%     plot(mean(rtmp,2),'DisplayName',num2str(lambdaList(i)));
%     hold on;
% end

% threshold is fixed to 0.05 to be conservative spike noise detection
% if strcmp(mode_name,'DoubleStep')
%     rr1 = logYifc_bprmvd - logYifc_model_bprmvd;
%     logYifc_bprmvd_isnan_new = or(logYifc_bprmvd_isnan,abs(rr1)>0.05);
%     outlier_added = abs(rr1)>0.05;
%     % remove nans
%     logYifc_bprmvd_new = interp_nan_column_given(logYifc_bprmvd,logYifc_bprmvd_isnan_new,logYifc_model_bprmvd);
% 
%     rr1_new = logYifc_bprmvd_new - logYifc_model_bprmvd;
%     % logYifc_bprmvd = interp_nan_column(logYifc_bprmvd,logYifc_bprmvd_isnan,wvc_bprmvd);
%     lambda_tmp = lambda_a * ( sum(abs(rr1_new(:))) ./ sum(abs(rr1(:))) ) ;
%     lambda_a_1(idxAlib) = lambda_tmp;
%     % perform again to remove
%     [ X1,Z1,C1,~,D1,rho,Rhov,~] = huwacbl1_gadmm_a(A_bprmvd,logYifc_bprmvd_new,wvc_bprmvd,'LAMBDA_A',lambda_a_1,...
%                                                    'CONCAVEBASE',C1,'Z0',Z1,'D0',D1,'X0',X1,'R0',rr1_new,'rho',rho,'Rhov',Rhov);
%     logBg_bprmvd = C1*Z1;
%     rr1 = logYifc_bprmvd - A_bprmvd*X1 - C1*Z1;
% 
%     outlier_added_new = abs(rr1)>0.05;
%     fprintf('%d,%d\n',sum(outlier_added(:)),sum(outlier_added_new(:)));
% end
% evaluate the leftover noise again

switch spike_detection_opt
%     case 'incremental'
%         rr1 = logYifc_bprmvd-logYifc_model_bprmvd;
%         logYifc_bprmvd_isnan = or(logYifc_bprmvd_isnan,abs(rr1)>0.05);
%         logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
    case 'original'
        rr1 = logYifc_bprmvd_ori-logYifc_model_bprmvd;
        rr1_std = nanstd(rr1,[],2);
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


R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;


Xlogtc_1d = sum(Xlogtc,1);
r_lad = zeros(Ny,L_bprmvd); d_lad = zeros(Ny,L_bprmvd); Rhov_lad = ones(Ny,1);
[logt_est_bprmvd,r_lad(vldpxl,:),d_lad(vldpxl,:),rho_lad,Rhov_lad(vldpxl,:)] = lad_gadmm_a_v2(Xlogtc_1d(:,vldpxl)', R_bprmvd(:,vldpxl)',...
                                'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad);
logt_est_bprmvd = logt_est_bprmvd';
%R_bprmvd_new = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
resNew_bprmvd = R_bprmvd - logt_est_bprmvd*Xlogtc_1d;
%resNewNrm_bprmvd = nansum(nansum(abs(resNew_bprmvd)));
resNewNrm_bprmvd = nanmedian(nanmedian(abs(resNew_bprmvd)));


%%
if isdebug
    spcs_bprmvd = logYifc_bprmvd - A_bprmvd(:,idxAlogtc)*Xlogtc;
    spcs_new = logYifc_bprmvd - logt_est_bprmvd*Xlogtc_1d;
    spc_r_bprmvd = logBg_bprmvd + logAB_bprmvd;
    logtc_best = nanmean((A_bprmvd(:,idxAlogtc)*Xlogtc)./sum(Xlogtc,1),2);
    diff_tList(:,1) = logt_est_bprmvd-logtc_best;
    logt_estList = [logtc_best logt_est_bprmvd];

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

    plot(ax_tr,wvc_bprmvd,logt_est_bprmvd,'DisplayName','iter=0');
    for k=84
        plot(ax_spc,wvc,exp(logYifc_cat(:,k)),'Color','k',...
            'DisplayName',sprintf('iter=0;%d cat\n',k));
        plot(ax_spc,wvc,exp(logYraifc_cat(:,k)),'Color',[0.5 0.5 0.5],...
            'DisplayName',sprintf('iter=0;%d cat\n',k));
        hold(ax_spc,'on');
        l1 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)),...
            'DisplayName',sprintf('iter=0;%d\n',k));
        l2 = plot(ax_spc,wvc_bprmvd,exp(spcs_new(:,k)),...
            'DisplayName',sprintf('iter=0;%d\n',k));
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
end

%% mainloop to improve the estimate
A_bprmvd = [logt_est_bprmvd Alib_bprmvd];
X = [Xlogtc_1d;Xlib];
% D = [sum(D1(idxAlogtc,:),1);D1(idxAlib,:)];
D = [zeros(1,size(D1,2)); D1(idxAlibstrt:end,:)];
C = C1;
Z = Z1;
lambda_a_2 = zeros([1,size(A_bprmvd,2)]);

rho = ones([1,Ny]);

% always update lambda_tmp
lambda_tmp = lambda_tmp*resNewNrm_bprmvd/resNrm_bprmvd;

for j=2:nIter+1
    if isdebug
        fprintf('Iter%d,lambda = %3.4e\n',j,lambda_tmp);
    end
    lambda_a_2(2:end) = lambda_tmp;   
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
    
    
    % detect spike noise
%     if strcmp(mode_name,'DoubleStep')
%         logYifc_model_bprmvd = A_bprmvd*X + C*Z;
%         rr = logYifc_bprmvd - logYifc_model_bprmvd;
%         logYifc_bprmvd_isnan_new = or(logYifc_bprmvd_isnan,abs(rr)>0.015);
%         % remove nans
%         logYifc_bprmvd_new = interp_nan_column_given(logYifc_bprmvd,logYifc_bprmvd_isnan_new,logYifc_model_bprmvd);
% 
%         rr_new = logYifc_bprmvd_new - logYifc_model_bprmvd;
%         lambda_tmp = lambda_tmp * ( sum(abs(rr_new(:))) ./ sum(abs(rr(:))) ) ;
%         lambda_a_2(idxAlibstrt:end) = lambda_tmp;
%         outlier_added = abs(rr)>0.015;
% 
%         [ X,Z,C,~,D,rho,Rhov ] = huwacbl1_gadmm_a(A_bprmvd,logYifc_bprmvd_new,wvc_bprmvd,...
%                                 'LAMBDA_A',lambda_a_2,'CONCAVEBASE',C,'Z0',Z,...
%                                 'D0',D,'X0',X,'R0',rr_new,...%'rho',rho,'Rhov',Rhov,...
%                                 'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb);
%         outlier_added_new = abs(rr)>0.015;
%         fprintf('%d,%d\n',sum(outlier_added(:)),sum(outlier_added_new(:)));
%     end
    
    logBg_bprmvd = C*Z;
    logAB_bprmvd = Alib_bprmvd*X(2:end,:);
    logYifc_model_bprmvd = logBg_bprmvd + logAB_bprmvd + A_bprmvd(:,1)*X(1,:);
    
    R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
    res_bprmvd = R_bprmvd-A_bprmvd(:,1)*X(1,:);
    % resNrm_bprmvd = nansum(nansum(abs(res_bprmvd)));
    resNrm_bprmvd = nanmedian(nanmedian(abs(res_bprmvd)));
    
    switch spike_detection_opt
%         case 'incremental'
%             rr = logYifc_bprmvd-logYifc_model_bprmvd;
%             logYifc_bprmvd_isnan = or(logYifc_bprmvd_isnan,abs(rr)>0.015);
%             logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
        case 'original'
            rr = logYifc_bprmvd_ori-logYifc_model_bprmvd;
            logYifc_bprmvd_isnan = or( logYifc_bprmvd_isnan_ori, abs(rr)>0.015 );
            logYifc_bprmvd = interp_nan_column_given(logYifc_bprmvd_ori,logYifc_bprmvd_isnan,logYifc_model_bprmvd);
            vldpxl = (sum(logYifc_bprmvd_isnan,1)/L_bprmvd) < 0.8;
        otherwise 
            error('not implemented');
    end   

    % remove nans
    
    
    %update logt_est!
    R_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
    [logt_est_bprmvd,r_lad(vldpxl,:),d_lad(vldpxl,:),rho_lad,Rhov_lad(vldpxl,:)]...
        = lad_gadmm_a_v2(X(1,vldpxl)', R_bprmvd(:,vldpxl)',...
        'X0',A_bprmvd(:,1)','R0',r_lad(vldpxl,:),'D0',d_lad(vldpxl,:),...
        'rho',rho_lad,'Rhov',Rhov_lad(vldpxl,:),...
        'tol',tol_lad,'maxiter',maxiter_lad,'verbose',verbose_lad);
    logt_est_bprmvd = logt_est_bprmvd';
%     logt_est = update_logt_est(R,X(1,:));
    diff_tList(:,j) = logt_est_bprmvd-logt_estList(:,j);
    
    if (sqrt(mean(diff_tList(:,j).^2,1))/sqrt(mean(diff_tList(:,1).^2,1)) < 0.05)
        % break;
    end
    
    
    %R_bprmvdNew = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd;
    resNew_bprmvd = R_bprmvd - logt_est_bprmvd*X(1,:);
    % resNewNrm_bprmvd = nansum(nansum(abs(resNew_bprmvd)));
    resNewNrm_bprmvd = nanmedian(nanmedian(abs(resNew_bprmvd)));
    
    A_bprmvd(:,1) = logt_est_bprmvd;
    
    lambda_tmp = lambda_tmp*resNewNrm_bprmvd/resNrm_bprmvd;

    if isdebug
        % figure(10); imsc(logYifc_bprmvd_isnan); drawnow;
    end
        
    if isdebug
        res_bprmvd = logYifc_bprmvd - logBg_bprmvd - logAB_bprmvd-A_bprmvd(:,1)*X(1,:);
        % resNrm_bprmvd = nansum(nansum(res_bprmvd.^2));
        RList(:,j) = vnorms(res_bprmvd,1,2);
        res_nrmList(j) = resNewNrm_bprmvd;
        logt_estList(:,j+1) = logt_est_bprmvd;
        diff_tList(:,j) = logt_est_bprmvd-logt_estList(:,j);
        spcs_bprmvd = logYifc_bprmvd - A_bprmvd(:,1)*X(1,:);
        spc_r_bprmvd = logBg_bprmvd + logAB_bprmvd;
        plot(ax_tr,wvc_bprmvd,logt_est_bprmvd,'DisplayName',sprintf('iter=%d',j));
        for k=84
            hold(ax_spc,'on');
            plot(ax_spc,wvc,exp(logYifc_cat(:,k)),'Color','k',...
                'DisplayName',sprintf('iter=0;%d cat\n',k));
            hold(ax_spc,'on');
            l1 = plot(ax_spc,wvc_bprmvd,exp(spcs_bprmvd(:,k)),...
                            'DisplayName',sprintf('iter=%d;%d\n',j,k));
            plot(ax_spc,wvc_bprmvd,exp(spc_r_bprmvd(:,k)),':','Color',l1.Color,...
                            'DisplayName',sprintf('iter=0;%d m\n',k));
            plot(ax_spc,wvc_bprmvd,exp(logBg_bprmvd(:,k)),'--','Color',l1.Color,...
                'DisplayName',sprintf('iter=%d;%d b\n',j,k));
        end
        hold(ax_spc,'off');
        plot(ax_dtr,wvc_bprmvd,diff_tList(:,j),'DisplayName',sprintf('iter=%d',j));
        plot(ax_res,res_nrmList);
        plot(ax_resv,wvc_bprmvd,res_bprmvd(:,:),'DisplayName',sprintf('iter=%d',j));
        plot(ax_mdtr,sqrt(mean(diff_tList.^2,1)));
        drawnow;
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

