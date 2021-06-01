function [Aicescalib,infoAicescalib,valid_idx] = crmsab_loadlibc_icescalib(opt,wabasename,c,bands_opt,wvc,varargin)

global localCRISM_PDSrootDir

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
Alibcachefname = constAicescalibcachefname(wabasename,opt,bands_opt,c,cntRmvl);
Alibcachefilepath = joinPath(Alibdir,Alibcachefname);
% Alibcachefilepath = joinPath(Alibdir,sprintf('Alib%s_%s_interp%s_b%s_c%03d.mat',libprefix,wabasename,optInterpid,bands_opt,c));
if ~overwrite && exist(Alibcachefilepath,'file')
    load(Alibcachefilepath,'Aicescalib');
    if nargout>1
        load(Alibcachefilepath,'infoAicescalib','valid_idx');
    end
else
    [Aicescalib,infoAicescalib] = wrapper_crism_icescalib(opt,wabasename,c);
    % select the band option
    bands = crmsab_genBands(bands_opt);
    Aicescalib= Aicescalib(bands,:) ; 
    % prun if there is any invalid one
    valid_idx = ~any(isnan(Aicescalib),1);
%     valid_idx2 = all(Aicelib>1e-10,1);
%     valid_idx = and(valid_idx1,valid_idx2);
    infoAicescalib = infoAicescalib(valid_idx);
    Aicescalib = Aicescalib(:,valid_idx);
    Aicescalib = normalizevec(Aicescalib,1);
    if cntRmvl
        [ continua,bases ] = learnContinuum(wvc,Aicescalib,'convhull');
        [ Aicescalib ] = CntRmvl(Aicescalib,'additive','CONTINUA',continua );
        Aicescalib = -Aicescalib;
    end
    
    save(Alibcachefilepath,'Aicescalib','infoAicescalib','bands','valid_idx');

end


end

function Alibcachefname = constAicescalibcachefname(wabasename,opt,bands_opt,c,cntRmvl)

Alibcachefname = sprintf('Aicelib%d_%s_b%d_c%03d_%d.mat',opt,wabasename,bands_opt,c,cntRmvl);

end