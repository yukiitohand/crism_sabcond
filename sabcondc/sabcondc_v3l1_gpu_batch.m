function [logYif_cor,logt_est,logAB,logBg,logIce,logYif_isnan,X] = sabcondc_v3l1_gpu_batch(logYif,WA,Alib,logT,BP,varargin)
% 
% INPUTS
%   logYif
%     [B x L x S] Each of the page is the observation
%   WA
%     [B x 1 x S]     Wavelength frame
%   Alib
%     [B x Nlib x S] Each of the page is library matrix for column s
%   logT
%     [B x Ntc x S] Each of the page is collection of transmission spectra
%     for column s
%   BP
%     [B x 1 x S]   Bad pixel information

% OUTPUTS

% [ logt_est,logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,X,vldpxl]

%%
%-------------------------------------------------------------------------%
% varargin
%-------------------------------------------------------------------------%
Aicelib = [];
maxiter_huwacb = int32(100);
tol_huwacb = 1e-4;
maxiter_lad = int32(1000);
tol_lad = 1e-4;
verbose_lad = 'no';
nIter = int32(5);
lambda_a = 0.01;
lambda_a_ice = 0;
verbose_huwacb = 'no';
precision = 'double';
th_badspc = 0.8;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'AICELIB'
                Aicelib = varargin{i+1};
            case 'NITER'
                nIter = varargin{i+1};
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
            case 'MAXITER_LAD'
                maxiter_lad = round(varargin{i+1});
                if (maxiter_huwacb <= 0 )
                       error('AL_iters must a positive integer');
                end
            case 'TOL_LAD'
                tol_lad = varargin{i+1};
            case 'VERBOSE_LAD'
                verbose_lad = varargin{i+1};
            case 'THRESHOLD_BADSPC'
                th_badspc = varargin{i+1};
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
[B,L,S]    = size(logYif);
[~,Nlib,~] = size(Alib);
[~,Nice,~] = size(Aicelib);
[~,Ntc,~]  = size(logT);

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
logYif  = gpuArray(logYif);
WA      = gpuArray(WA);
Alib    = gpuArray(Alib);
Aicelib = gpuArray(Aicelib);
logT    = gpuArray(logT);
BP      = gpuArray(BP);
logYif_isnan = gpuArray(logYif_isnan);
logYif_isnan_ori = gpuArray(logYif_isnan_ori);
badspc = gpuArray(badspc);

%%
%-------------------------------------------------------------------------%
% initialization using the matrix logT
%-------------------------------------------------------------------------%
A = cat(2,logT,Alib);
% compute concave basese beforehand
% C = concaveOperator(WA);
% Dinv = C \ eye(B);
C = concaveOperator_v2(WA);
Dinv = pagefun(@mldivide,C,eye(B,'gpuArray'));
s_d = vnorms(Dinv,1);
C = Dinv./s_d;
clear Dinv s_d;
[~,Nc,~] = size(C);
if strcmpi(precision,'single')
    C = single(C);
end

% create lambdas
% lambda_c = zeros(B,L,S,precision);
% lambda_a2 = zeros(Nlib+Ntc,L,S,precision);
% lambda_r = ones(B,L,S,precision);
lambda_c = zeros(B,L,S,precision,'gpuArray');
lambda_a_2 = zeros(Nlib+Nice+Ntc,L,S,precision,'gpuArray');
lambda_r = ones(B,L,S,precision,'gpuArray');

lambda_a_2(:,(1+Ntc):(Ntc+Nice)) = lambda_a_ice;
lambda_a_2(:,(1+Ntc+Nice):end,:) = lambda_a;

% 
lambda_c(logYif_isnan) = inf;
lambda_c([1,Nc],:,:) = 0; % safeguard
lambda_r(logYif_isnan) = 0;

c2_z = zeros([Nc,1],precision);
% c2_z = zeros([Nc,1],precision,'gpuArray');
c2_z(1) = -inf; c2_z(Nc) = -inf;

% main computation
% tic; 
[ X,Z,~,D,rho,Rhov,~,~,cost_val]...
    = huwacbl1_admm_gat_a_batch(A,logYif,C,'LAMBDA_A',lambda_a_2,...
            'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,'C2_Z',c2_z,......
            'tol',1e-5,'maxiter',maxiter_huwacb,'verbose',verbose_huwacb,...
            'precision',precision); % toc;
% tic;
% [ X1,Z1,~,R,D1,rho1,Rhov1,~,~,cost_val1]...
% = huwacbl1_admm_gat_a(A,logYif,gather(WA),'LAMBDA_A',lambda_a2,...
% 'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
% 'tol',1e-5,'maxiter',100,'verbose','no',...
% 'precision',precision,'gpu',true,'Concavebase',C,'debug',false);
% toc;

% evaluate bad pixels
RR = logYif - pagefun(@mtimes,A,X) - pagefun(@mtimes,C,Z);
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
RR  = RR + pagefun(@mtimes,logT,X(1:Ntc,:,:));
Xtc = sum(X(1:Ntc,:,:),1);
[logt_est,r_lad,d_lad,rho_lad,Rhov_lad,~,~,cost_val,Kcond]...
       = lad_admm_gat_b_batch(permute(Xtc,[2,1,3]), permute(RR,[2,1,3]),...
            'lambda_r',permute(lambda_r,[2,1,3]),'tol',tol_lad,'maxiter',maxiter_lad,...
            'verbose',verbose_lad,'precision',precision);
%
logt_est = permute(logt_est,[2,1,3]);
RR = RR - pagefun(@mtimes,logt_est,Xtc);
resNewNrm = nansum(abs(lambda_r .* RR),[1,2]);

%%
%-------------------------------------------------------------------------%
% main loop
%-------------------------------------------------------------------------%
A = cat(2,logt_est,Alib);
X = cat(1,Xtc,X(1+Ntc:end,:,:));
D = cat(1,zeros(1,L,S,precision,'gpuArray'),D(1+Ntc:end,:,:));
lambda_a_2 = ones(1+Nice+Nlib,L,S,precision,'gpuArray');
lambda_a_2(1,:,:) = 0;
lambda_a_2(:,2:(Nice+1)) = lambda_a_ice;
lambda_a_2(:,(2+Nice):end,:) = lambda_a;
% rho = ones([1,L,S],precision,'gpuArray');
Rhov = cat(1,ones(1,1,S,precision,'gpuArray'),Rhov(Ntc+1:Ntc+Nlib+Nc+B,:,:));
% lambda_tmp = lambda_a;
% always update lambda_tmp
% lambda_tmp = lambda_tmp .* resNewNrm ./ resNrm;

for n=2:nIter
    lambda_a_2(2:end,:,:) = lambda_a_2(2:end,:,:) .* resNewNrm ./ resNrm;
    % tic;
    if n==2
        [ X,Z,~,D,rho,Rhov,~,~,cost_val,Tcond ]...
            = huwacbl1_admm_gat_a_batch(A,logYif,C,...
                'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
                'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
                'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb,...
                'PRECISION',precision);
%         [ X1,Z1,~,R,D1,rho1,Rhov1,~,~,cost_val1]...
%         = huwacbl1_admm_gat_a(A,logYif,gather(WA),'LAMBDA_A',lambda_a2,...
%         'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,.....
%         'Z0',Z,'D0',D,'X0',X,'R0',RR,...
%         'tol',1e-5,'maxiter',100,'verbose','no',...
%         'precision',precision,'gpu',true,'Concavebase',C,'debug',false);
    else
       [ X,Z,~,D,rho,Rhov,~,~,cost_val ] = huwacbl1_admm_gat_a_batch(A,logYif,C,...
          'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
          'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
          'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
          'PRECISION',precision,'Tcond',Tcond);
    end
    % toc;
    % evaluate bad pixels
    RR = logYif - pagefun(@mtimes,A,X) - pagefun(@mtimes,C,Z);
    logYif_isnan = or( logYif_isnan_ori, abs(RR)>0.015 );
    badspc = (sum(logYif_isnan,1)/B) > th_badspc;
    logYif_isnan = or(logYif_isnan,badspc);

    lambda_c(logYif_isnan) = inf; lambda_c(~logYif_isnan) = 0;
    lambda_c([1,Nc],:,:) = 0; % safeguard
    lambda_r(logYif_isnan) = 0; lambda_r(~logYif_isnan) = 1;
    
    resNrm = nansum(abs(lambda_r .* RR),[1,2]);
    
    % update logt_est
    RR  = RR + pagefun(@mtimes,A(:,1,:),X(1,:,:));
    [logt_est,r_lad,d_lad,rho_lad,Rhov_lad,~,~,cost_val]...
       = lad_admm_gat_b_batch(permute(Xtc,[2,1,3]), permute(RR,[2,1,3]),...% 'rho',rho_lad,'Rhov',Rhov_lad,...
            'lambda_r',permute(lambda_r,[2,1,3]),'tol',tol_lad,'maxiter',maxiter_lad,...
            'verbose',verbose_lad,'precision',precision,'Kcond',Kcond);
    logt_est = permute(logt_est,[2,1,3]);
    RR = RR - pagefun(@mtimes,logt_est,X(1,:,:));
    resNewNrm = nansum(abs(lambda_r .* RR),[1,2]);
    
    A(:,1,:) = logt_est;
    
    % lambda_tmp = lambda_tmp*resNewNrm/resNrm;
    
end

%%
%-------------------------------------------------------------------------%
% last iteration
%-------------------------------------------------------------------------%
[ X,Z,D,rho,Rhov,~,~,cost_val ] = huwacbl1_admm_gat_a_batch(A,logYif,C,...
      'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
      'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,...%'rho',rho,'Rhov',Rhov,...
      'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
      'PRECISION',precision);

% substituting all the variables
logAB        = pagefun(@mtimes,Alib,X((2+Nice):end,:,:));
logBg        = pagefun(@mtimes,C,Z);
logIce       = pagefun(@mtimes,Aicelib,X(2:(Nice+1),:,:));
logt_est     = A(:,1,:);
logYif_cor   = logYif - pagefun(@mtimes,logt_est,X(1,:,:)) - logIce;
logYif_cor(logYif_isnan) = nan;

[logYif_cor,logt_est,logAB,logBg,logYif_isnan,X] = gather(logYif_cor,logt_est,logAB,logBg,logYif_isnan,X);


end

