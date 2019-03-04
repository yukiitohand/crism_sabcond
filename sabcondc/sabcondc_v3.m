function [ logt_est,logYifc_cor,logAB,logBg,ancillary] = sabcondc_v3( Alib,logYifc,wvc,logtc,varargin )
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

maxiter_huwacb = 3000;
tol_huwacb = 1e-5;
nIter = 0;
vis = 0;
lambda_a = 0.001;
logYifc_cat = 0;
logYraifc_cat = 0;
verbose_huwacb = 'no';
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
            case 'VERBOSE'
                verbose = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

A = [logtc Alib];
idxAlibstrt = size(logtc,2)+1;
lambda_a_1 = zeros([1,size(A,2)]);
lambda_a_1(idxAlibstrt:end) = lambda_a;

% 
% logYifc(126:128,32:33) = nan;
% logYifc(126:128,190:192) = nan;
logYifc_ori = logYifc;
logYifc_bprmvd = logYifc .* gp;
[ groups ] = nan_grouper( logYifc_bprmvd );

%% initialization of logt_est
[ X1,logBg,~,~,ancillaries ] = nanhuwacb(A,logYifc_bprmvd,wvc,groups,...
                'MAXITER',maxiter_huwacb,'TOL',tol_huwacb,'VERBOSE',verbose_huwacb,...
                'LAMBDA_A',lambda_a_1);

X = X1(idxAlibstrt:end,:); Xlogtc = X1(1:idxAlibstrt-1,:);
logAB = Alib*X;
R = logYifc - logBg - logAB;
res = R - logtc * Xlogtc;
resNrm = nansum(nansum(res.^2));

logt_est = R ./ sum(Xlogtc,1);
logt_est = nanmean(logt_est,2);
  
resNew = R - logt_est*sum(Xlogtc,1);
resNewNrm = nansum(nansum(resNew.^2));

%%
if vis
    spcs = logYifc - A(:,1:idxAlibstrt-1)*X1(1:idxAlibstrt-1,:);
    spcs_bprmvd = logYifc_bprmvd - A(:,1:idxAlibstrt-1)*X1(1:idxAlibstrt-1,:);
    spc_r = logBg + logAB;
    logtc_best = nanmean((A(:,1:idxAlibstrt-1)*X1(1:idxAlibstrt-1,:))./sum(X1(1:(idxAlibstrt-1),:),1),2);
    diff_tList(:,1) = logt_est-logtc_best;
    logt_estList = [logtc_best logt_est];

    RList(:,1) = vnorms(res,1,2);
    res_nrmList(1) = resNrm;

    figure; ax_tr = subplot(1,1,1);movegui(gcf,'northwest');
    plot(ax_tr,wvc,A(:,1:idxAlibstrt-1),'Color',[0.5 0.5 0.5],'DisplayName','iter=0'); hold(ax_tr,'on'); 
    figure; ax_dtr = subplot(1,1,1); hold(ax_dtr,'on');movegui(gcf,'north');
    figure; ax_spc = subplot(1,1,1); hold(ax_spc,'off');movegui(gcf,'northeast');
    figure; ax_mdtr = subplot(1,1,1); movegui(gcf,'southwest');
    title(ax_mdtr,'norm of difference of consequtive logt');
    figure; ax_res = subplot(1,1,1); movegui(gcf,'south');
    title(ax_mdtr,'norm of residual');
    figure; ax_resv = subplot(1,1,1); hold(ax_resv,'off');movegui(gcf,'southeast');

    plot(ax_tr,wvc,logt_est,'DisplayName','iter=0');
    for k=300
        plot(ax_spc,wvc,exp(logYifc_cat(:,k)),'Color','k',...
            'DisplayName',sprintf('iter=0;%d cat\n',k));
        hold(ax_spc,'on');
        l1 = plot(ax_spc,wvc,exp(spcs(:,k)),...
            'DisplayName',sprintf('iter=0;%d\n',k));

%         plot(ax_spc,wvc,exp(spc_r(:,k)),'-','Color',l1.Color,...
%             'DisplayName',sprintf('iter=0;%d m\n',k));
%         plot(ax_spc,wvc,exp(logBg(:,k)),'Color',l1.Color,...
%             'DisplayName',sprintf('iter=0;%d b\n',k));
        hold off
%             pause

    end
    hold(ax_spc,'off');
    plot(ax_dtr,wvc,diff_tList(:,1),'DisplayName','iter=0');
    plot(ax_res,res_nrmList);
    plot(ax_resv,wvc,res(:,:),'DisplayName','iter=0');
    plot(ax_mdtr,sqrt(mean(diff_tList.^2,1)));

    drawnow;
end

%% mainloop to improve the estimate
A = [logt_est Alib];
X = [nansum(Xlogtc,1);X];
lambda_a_2 = zeros([1,size(A,2)]);

if resNewNrm/resNrm<0.5 
    %update the trade-off parameter based on improvement in fitting
    lambda_tmp = lambda_a*resNewNrm/resNrm;
else
    lambda_tmp = lambda_a;
end

for j=2:nIter+1
%     fprintf('Iter%d,lambda = %3.4e\n',j,lambda_tmp);
    lambda_a_2(2:end) = lambda_tmp;   
    logt_estList(:,j) = A(:,1);
    
    for gidi = 1:length(groups)
        activeIdxes = groups(gidi).idxs;
        activeBands = ~groups(gidi).ptrn;
        if mean(activeBands)>0.3 % only perform when more than 30% are good bands.
            ytmp = logYifc(activeBands,activeIdxes);
            Atmp = A(activeBands,:);
            wvtmp = wvc(activeBands);
            if j==2
                Dtmp = ancillaries(gidi).D;
                Dtmp = [sum(Dtmp(1:idxAlibstrt-1,:),1);Dtmp(idxAlibstrt:end,:)];
            else
                Dtmp = ancillaries(gidi).D;
            end
            
            [xtmp,ztmp,Ctmp,dtmp,rho] = huwacb_admm2_test(Atmp,ytmp,wvtmp,...
                'TOL',tol_huwacb,'MAXITER',maxiter_huwacb,'LAMBDA_A',lambda_a_2,...
                'CONCAVEBASE',ancillaries(gidi).C,'Z0',ancillaries(gidi).Z,...
                'D0',Dtmp,'X0',X(:,activeIdxes),'rho',ancillaries(gidi).rho,...
                'verbose',verbose_huwacb);

            X(:,activeIdxes) = xtmp;
            logBgtmp = Ctmp * ztmp;
            logABtmp = Atmp(:,2:end) * xtmp(2:end,:);

            logBg(activeBands,activeIdxes) = logBgtmp;
            logAB(activeBands,activeIdxes) = logABtmp;
            
            ancillaries(gidi).activeIdxes = activeIdxes;
            ancillaries(gidi).activeBands = activeBands;
            ancillaries(gidi).C = Ctmp;
            ancillaries(gidi).Z = ztmp;
            ancillaries(gidi).rho = rho;
            ancillaries(gidi).D = dtmp;

        end
    end

    %update logt_est. nan values are ignored.
    R = logYifc - logBg - logAB;
    logt_est = update_logt_est(R,X(1,:));
    diff_tList(:,j) = logt_est-logt_estList(:,j);
    
    if (sqrt(mean(diff_tList(:,j).^2,1))/sqrt(mean(diff_tList(:,1).^2,1)) < 0.05)
%         break;
    end
    
    res = R-A(:,1)*X(1,:);
    resNrm = nansum(nansum(res.^2));
    resNew = R - logt_est*X(1,:);
    resNewNrm = nansum(nansum(resNew.^2));
    
    A(:,1) = logt_est;
    
    if resNewNrm/resNrm<0.5
        lambda_tmp = lambda_tmp*resNewNrm/resNrm;
    end
        
    if vis
        res = logYifc - logBg - logAB-A(:,1)*X(1,:);
        resNrm = nansum(nansum(res.^2));
        RList(:,j) = vnorms(res,1,2);
        res_nrmList(j) = norm(resNrm,'fro');
        logt_estList(:,j+1) = logt_est;
        diff_tList(:,j) = logt_est-logt_estList(:,j);
        spcs = logYifc - A(:,1)*X(1,:);
        spc_r = logBg + logAB;
        plot(ax_tr,wvc,logt_est,'DisplayName',sprintf('iter=%d',j));
        for k=300
            hold(ax_spc,'on');
            plot(ax_spc,wvc,exp(logYifc_cat(:,k)),'Color','k',...
                'DisplayName',sprintf('iter=0;%d cat\n',k));
            hold(ax_spc,'on');
            l1 = plot(ax_spc,wvc,exp(spcs(:,k)),...
                            'DisplayName',sprintf('iter=%d;%d\n',j,k));
%             plot(ax_spc,wvc,exp(spc_r(:,k)),'Color',l1.Color,...
%                             'DisplayName',sprintf('iter=0;%d m\n',k));
%             plot(ax_spc,wvc,exp(logBg(:,k)),'Color',l1.Color,...
%                 'DisplayName',sprintf('iter=%d;%d b\n',j,k));
            hold(ax_spc,'off');
        end
        hold(ax_spc,'off');
        plot(ax_dtr,wvc,diff_tList(:,j),'DisplayName',sprintf('iter=%d',j));
        plot(ax_res,res_nrmList);
        plot(ax_resv,wvc,res(:,:),'DisplayName',sprintf('iter=%d',j));
        plot(ax_mdtr,sqrt(mean(diff_tList.^2,1)));
        drawnow;
    end
end

%% last iteration after estimating log_est
for gidi = 1:length(groups)
    activeIdxes = groups(gidi).idxs;
    activeBands = ~groups(gidi).ptrn;
    if mean(activeBands)>0.3 % only perform when more than 30% are good bands.
        ytmp = logYifc(activeBands,activeIdxes);
        Atmp = A(activeBands,:);
        wvtmp = wvc(activeBands);

        Dtmp = ancillaries(gidi).D;

        [xtmp,ztmp,Ctmp,dtmp,rho] = huwacb_admm2_test(Atmp,ytmp,wvtmp,...
            'CONCAVEBASE',ancillaries(gidi).C,'Z0',ancillaries(gidi).Z,...
            'D0',Dtmp,'X0',X(:,activeIdxes),'rho',ancillaries(gidi).rho,...
            'TOL',tol_huwacb,'MAXITER',maxiter_huwacb,'LAMBDA_A',lambda_a_2,...
            'verbose',verbose_huwacb);
        % maybe lambda_a_2 instead of lambda_a is right(2017/10/11, Yuki)
%         [xtmp,ztmp,Ctmp,dtmp,rho] = huwacb_admm2_test(Atmp,ytmp,wvtmp,...
%             'CONCAVEBASE',ancillaries(gidi).C,'Z0',ancillaries(gidi).Z,...
%             'D0',Dtmp,'X0',X(:,activeIdxes),'rho',ancillaries(gidi).rho,...
%             'TOL',tol_huwacb,'MAXITER',maxiter_huwacb,'LAMBDA_A',lambda_a,...
%             'verbose',verbose_huwacb);

        logBgtmp = Ctmp * ztmp;
        logABtmp = Atmp(:,2:end) * xtmp(2:end,:);

        X(:,activeIdxes) = xtmp;

        logBg(activeBands,activeIdxes) = logBgtmp;
        logAB(activeBands,activeIdxes) = logABtmp;

        ancillaries(gidi).activeIdxes = activeIdxes;
        ancillaries(gidi).activeBands = activeBands;
        ancillaries(gidi).C = Ctmp;
        ancillaries(gidi).Z = ztmp;
        ancillaries(gidi).rho = rho;
        ancillaries(gidi).D = dtmp;

    end
end

logYifc_cor = logYifc - logt_est*X(1,:);

ancillary = [];
ancillary.X = X;
ancillary.lambda.init = lambda_a; % edited by Yuki 2017/10/11
ancillary.lambda.last = lambda_a_2(2:end); % edited by Yuki 2017/10/11, the final lambda used
ancillary.nIter = nIter;
ancillary.huwacb_func = 'huwacb_admm2';
ancillary.maxiter_huwacb = maxiter_huwacb;
ancillary.tol_huwacb = tol_huwacb;
% ancillary.opt_tinit = opt_tinit;


end

