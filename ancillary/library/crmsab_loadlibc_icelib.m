function [Aicelib,infoAicelib,valid_idx] = loadlibc_crism_icelib(opt,wabasename,c,bands_opt,wvc,varargin)

global crism_env_vars
localCRISM_PDSrootDir = crism_env_vars.localCRISM_PDSrootDir;

overwrite = 0;
cntRmvl = 0;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'CNTRMVL'
                cntRmvl = varargin{i+1};
            case 'OVERWRITE'
                overwrite = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

Alibdir = joinPath(localCRISM_PDSrootDir,'cache/WA/',wabasename);
Alibcachefname = constAicelibcachefname(wabasename,opt,bands_opt,c,cntRmvl);
Alibcachefilepath = joinPath(Alibdir,Alibcachefname);
% Alibcachefilepath = joinPath(Alibdir,sprintf('Alib%s_%s_interp%s_b%s_c%03d.mat',libprefix,wabasename,optInterpid,bands_opt,c));
if ~overwrite && exist(Alibcachefilepath,'file')
    load(Alibcachefilepath,'Aicelib');
    if nargout>1
        load(Alibcachefilepath,'infoAicelib','valid_idx');
    end
else
    [Aicelib,infoAicelib] = wrapper_crism_icelib(opt,wabasename,c);
    % select the band option
    bands = genBands(bands_opt);
    Aicelib= Aicelib(bands,:) ; 
    % prun if there is any invalid one
    valid_idx = ~any(isnan(Aicelib),1);
%     valid_idx2 = all(Aicelib>1e-10,1);
%     valid_idx = and(valid_idx1,valid_idx2);
    infoAicelib = infoAicelib(valid_idx);
    Aicelib = Aicelib(:,valid_idx);
    if cntRmvl
        [ continua,bases ] = learnContinuum(wvc,Aicelib,'convhull');
        [ Aicelib ] = CntRmvl(Aicelib,'additive','CONTINUA',continua );
        Aicelib = -Aicelib;
    end
    
    save(Alibcachefilepath,'Aicelib','infoAicelib','bands','valid_idx');

end


end

