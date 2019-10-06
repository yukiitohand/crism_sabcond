function [logYif_cor,logt_est,logAB,logBg,logYif_isnan,X] = sabcondc_v3l1_gpu_batch(logYif,WA,Alib,logT,BP,varargin)
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

logYif = gpuArray(logYif);
WA     = gpuArray(WA);
Alib   = gpuArray(Alib);
logT   = gpuArray(logT);
BP     = gpuArray(BP);

% [ logt_est,logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,X,vldpxl]

%%
%-------------------------------------------------------------------------%
% varargin
%-------------------------------------------------------------------------%
maxiter_huwacb = int32(100);
tol_huwacb = 1e-4;
maxiter_lad = int32(1000);
tol_lad = 1e-4;
verbose_lad = 'no';
nIter = int32(5);
lambda_a = 0.01;
verbose_huwacb = 'no';
precision = 'double';
th_badspc = 0.8;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
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
            case 'THRESHOLD_BADSPC'
                th_badspc = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end




%%
%-------------------------------------------------------------------------%
% Preprocessing
%-------------------------------------------------------------------------%
[B,L,S]    = size(logYif);
[~,Nlib,~] = size(Alib);
[~,Ntc,~]  = size(logT);

% preprocessing is necessary because some of the values in logYif might be
% nans due to harse noise.
logYif_isnan = isnan(logYif);
% replace those values with linear interpolation
logYif = interp_nan_column(logYif,logYif_isnan,WA);
% mark bad spectra
badspc = (sum(logYif_isnan,1)/B) > 0.5;

% mark bad pixels and badspc
logYif_isnan = or(logYif_isnan,BP);
logYif_isnan = or(logYif_isnan,badspc);

% save the original bad entries for future use.
logYif_isnan_ori = logYif_isnan;
for i=1:S
    logYif(:,:,i) = interp_nan_column(logYif(:,:,i),logYif_isnan(:,:,i),WA(:,i));
end

%%
%-------------------------------------------------------------------------%
% initialization using the matrix logT
%-------------------------------------------------------------------------%
A = cat(2,logT,Alib);
% compute concave basese beforehand
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
lambda_c = zeros(B,L,S,precision,'gpuArray');
lambda_a2 = zeros(Nlib+Ntc,L,S,precision,'gpuArray');
lambda_r = ones(B,L,S,precision,'gpuArray');

% 
lambda_c(logYif_isnan) = inf;
lambda_r(logYif_isnan) = 0;

c2_z = zeros([Nc,1],precision,'gpuArray');
c2_z(1) = -inf; c2_z(Nc) = -inf;

% main computation
[ X,Z,D,rho,Rhov,~,~,cost_val]...
    = huwacbl1_admm_gat_a_batch(A,logYif,C,'LAMBDA_A',lambda_a2,...
            'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,'C2_Z',c2_z,......
            'tol',1e-5,'maxiter',maxiter_huwacb,'verbose',verbose_huwacb,...
            'precision',precision);

% evaluate bad pixels
RR = logYif - pagefun(@mtimes,A,X) - pagefun(@mtimes,C,Z);
RR_std = nanstd(RR,[],2);
logYif_isnan = or( logYif_isnan_ori, abs(RR)>0.1 );
logYif_isnan = or(logYif_isnan,RR_std>0.015);
% finally flag spectra that have too many bad channels.
badspc = (sum(logYif_isnan,1)/B) > th_badspc;
logYif_isnan = or(logYif_isnan,badspc);


lambda_c(logYif_isnan) = inf; lambda_c(~logYif_isnan) = 0;
lambda_r(logYif_isnan) = 0; lambda_r(~logYif_isnan) = 1;

resNrm = nansum(abs(RR).* lambda_r,'all');

% Get initial transmission spectrum
RR  = RR + pagefun(@mtimes,logT,X(1:Ntc,:,:));
Xtc = sum(X(1:Ntc,:,:),1);
[logt_est,r_lad,d_lad,rho_lad,Rhov_lad,~,~,cost_val]...
       = lad_admm_gat_b_batch(permute(Xtc,[2,1,3]), permute(RR,[2,1,3]),...
            'lambda_r',lambda_r','tol',tol_lad,'maxiter',maxiter_lad,...
            'verbose',verbose_lad,'precision',precision);
%
logt_est = permute(logt_est,[2,1,3]);
RR = RR - pagefun(@mtimes,logt_est,Xtc);
resNewNrm = nansum(abs(lambda_r .* RR),'all');

%%
%-------------------------------------------------------------------------%
% main loop
%-------------------------------------------------------------------------%
A = cat(2,logt_est,Alib);
X = cat(1,Xtc,X(idxAlib,:,:));
D = cat(1,zeros(1,L,S,precision,'gpuArray'),D(idxAlibstrt:end,:,:));
lambda_a_2 = zeros(1,1+Nlib,M,precision,'gpuArray');
rho = ones([1,Ny,M],precision,'gpuArray');
lambda_tmp = lambda_a;
% always update lambda_tmp
lambda_tmp = lambda_tmp*resNewNrm/resNrm;

for n=1:nIter
    lambda_a_2(2:end) = lambda_tmp;
    if n==2
        [ X,Z,C,~,D,rho,Rhov ]...
            = huwacbl1_admm_gat_(A,logYif,C,...
                'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
                'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,...
                'verbose',verbose_huwacb,'tol',1e-5,'maxiter',maxiter_huwacb,...
                'PRECISION',precision);
    else
       [ X,Z,C,~,D,rho,Rhov ] = huwacbl1_admm_gat_a(A,logYif,C,...
          'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
          'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
          'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
          'PRECISION',precision);
    end
    
    % evaluate bad pixels
    RR = logYif - pagefun(@mtimes,A,X) - pagefun(@mtimes,C,Z);
    logYif_isnan = or( logYif_isnan_ori, abs(RR)>0.015 );
    badspc = (sum(logYif_isnan,1)/B) > th_badspc;
    logYif_isnan = or(logYif_isnan,badspc);

    lambda_c(logYif_isnan) = inf; lambda_c(~logYif_isnan) = 0;
    lambda_r(logYif_isnan) = 0; lambda_r(~logYif_isnan) = 1;
    
    resNrm = nansum(abs(lambda_r .* RR),'all');
    
    % update logt_est
    RR  = RR + pagefun(@mtimes,A(:,1,:),X(1,:,:));
    [logt_est,r_lad,d_lad,rho_lad,Rhov_lad,~,~,cost_val]...
       = lad_admm_gat_b_batch(permute(Xtc,[2,1,3]), permute(RR,[2,1,3]),...
            'X0',A(:,1)','R0',r_lad,'D0',d_lad,'rho',rho_lad,'Rhov',Rhov_lad,...
            'lambda_r',lambda_r,'tol',tol_lad,'maxiter',maxiter_lad,...
            'verbose',verbose_lad,'precision',precision);
    logt_est = permute(logt_est,[2,1,3]);
    RR = RR - pagefun(@mtimes,logt_est,X(1,:,:));
    resNewNrm = nansum(abs(lambda_r .* RR),'all');
    
    A(:,1,:) = logt_est;
    
    lambda_tmp = lambda_tmp*resNewNrm/resNrm;
    
end

%%
%-------------------------------------------------------------------------%
% last iteration
%-------------------------------------------------------------------------%
[ X,Z,C,~,D,rho,Rhov ] = huwacbl1_admm_gat_a(A,logYif,C,...
      'LAMBDA_A',lambda_a_2,'LAMBDA_C',lambda_c,'LAMBDA_R',lambda_r,...
      'C2_Z',c2_z,'Z0',Z,'D0',D,'X0',X,'R0',RR,'rho',rho,'Rhov',Rhov,...
      'verbose',verbose_huwacb,'tol',tol_huwacb,'maxiter',maxiter_huwacb,...
      'PRECISION',precision);

% substituting all the variables
logAB        = pagefun(@mtimes,Alib,X(2:end,:,:));
logBg        = pagefun(@mtimes,C,Z);
logt_est     = A(:,1,:);
logYif_cor   = logYif - pagefun(logt_est,X(1,:,:));


end

