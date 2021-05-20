function [Ahitranlib,infoAhitranlib,valid_idx] = loadlibc_crism_hitranlib(opt,wabasename,c,bands_opt,wvc,varargin)

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
Alibcachefname = constAhitranlibcachefname(wabasename,opt,bands_opt,c,cntRmvl);
Alibcachefilepath = joinPath(Alibdir,Alibcachefname);
% Alibcachefilepath = joinPath(Alibdir,sprintf('Alib%s_%s_interp%s_b%s_c%03d.mat',libprefix,wabasename,optInterpid,bands_opt,c));
if ~overwrite && exist(Alibcachefilepath,'file')
    load(Alibcachefilepath,'Ahitranlib');
    if nargout>1
        load(Alibcachefilepath,'infoAhitranlib','valid_idx');
    end
else
    [Ahitranlib,infoAhitranlib] = wrapper_crism_hitranlib(opt,wabasename,c);
    % select the band option
    crmsab_bands = crmsab_genBands(bands_opt);
    Ahitranlib = Ahitranlib(bands,:); 
    % prun if there is any invalid one
    valid_idx = ~any(isnan(Ahitranlib),1);
%     valid_idx2 = all(Ahitranlib>1e-10,1);
%     valid_idx = and(valid_idx1,valid_idx2);
    infoAhitranlib = infoAhitranlib(valid_idx);
    Ahitranlib = Ahitranlib(:,valid_idx);
    if cntRmvl
        [ continua,bases ] = learnContinuum(wvc,Ahitranlib,'convhull');
        [ Ahitranlib ] = CntRmvl(Ahitranlib,'additive','CONTINUA',continua );
        Ahitranlib = -Ahitranlib;
    end
    
    save(Alibcachefilepath,'Ahitranlib','infoAhitranlib','bands','valid_idx');

end


end

function Alibcachefname = constAhitranlibcachefname(wabasename,opt,bands_opt,c,cntRmvl)

Alibcachefname = sprintf('Ahitranlib%d_%s_b%d_c%03d_%d.mat',opt,wabasename,bands_opt,c,cntRmvl);

end