function [ info_sabcond ] = const_sabcond_info(TRRIFdata,varargin )
% [ infohsi ] = loadinfo_procHSI_v3( crsc,varargin )
% load information of data processed with sabcond
%   Input Parameters
%      TRRIFdata: CRISM data used for processing 
%   Optional Parameters
%      'METHODTYPE' : method type for the product
%                     (default) 'sabcondv2'
%      'LIBPREFIX': prefix for the library
%                   (default) 'Lib1112'
%      'SAVE_DIR':  directory for the image cube
%                   (default) './resu/cccnnnnnn/'
%                            ccc = class type of the obervation
%                            nnnnnnn = observation id
%      'BANDS_OPT': band information used for the processing,
%                   option index.
%                   (default) 3
%      'OPTINTERPID': option id for how the library is interpolated
%                     (default) 1
%      'NITER':     number of outerloop iterations
%                   (default) 20
%      'CNTRMVL': continuum removal option
%                  (default) 1
%      'CNTRMVL_ICE': continuum removal option for ice lib
%                  (default) 1
%      'ADDITIONAL_SUFFIX': suffix at the end of the name
%                           (default) ''
%      'ALIBDIR': dir path to the cached library
%                 (default)[localCRISM_PDSrootDir]/cache/WA/[Basename of WA]/
%                 
%   Output Parameters
%      info_sabcond: struct, having the fields:
%         Alibdir
%         fname_supple
%         libprefix
%         optInterpid
%         bands_opt
%         cntRmvl
%         basenameWA
%         wa
%         bands
%         lines
%         T_est;
%         basename
%         save_dir
%
%   Notes:

global crism_env_vars
localCRISM_PDSrootDir = crism_env_vars.localCRISM_PDSrootDir;
mt = 'sabcondv2la2';
libprefix = 'Lib1112';
save_dir = sprintf('./resu/%s/',TRRIFdata.dirname);
bands_opt = 3;
optInterpid = 1;
nIter = 20;
cntRmvl = 1;
cntRmvl_ice = 0;
additional_suffix = '';

if ~isfield(TRRIFdata.basenamesCDR,'WA')
    TRRIFdata.load_basenamesCDR();
end

Alibdir = joinPath(localCRISM_PDSrootDir,'cache/WA/',TRRIFdata.basenamesCDR.WA);

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'METHODTYPE'
                mt = varargin{i+1};
            case 'LIBPREFIX'
                libprefix = varargin{i+1};
            case 'SAVE_DIR'
                save_dir = varargin{i+1};
            case 'BANDS_OPT'
                bands_opt = varargin{i+1};
            case 'OPTINTERPID'
                setups = [];
                optInterpid = varargin{i+1};
            case 'NITER'
                nIter = varargin{i+1};
            case 'CNTRMVL'
                cntRmvl = varargin{i+1};
            case 'CNTRMVL_ICE'
                cntRmvl_ice = varargin{i+1};
            case 'ADDITIONAL_SUFFIX'
                additional_suffix = varargin{i+1};
            case 'ALIBDIR'
                Alibdir = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

[suffix] = const_suffix_v2(mt,cntRmvl,libprefix,optInterpid,bands_opt,nIter,additional_suffix);
[optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib] = decompose_libprefix(libprefix);
optLibs = [optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib];
libprefix_Alib = const_libprefix(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib);
if isempty(opticelib)
    cntRmvl_ice = '';
end

%
basename_cr = [TRRIFdata.basename suffix];
fname_supple = joinPath(save_dir, [TRRIFdata.basename suffix '.mat']);
load(fname_supple,'wa','bands','T_est');

% wa = crim.img.cdr.WA;
% wa = squeeze(wa)';

info_sabcond = [];
info_sabcond.Alibdir = Alibdir;
info_sabcond.fname_supple = fname_supple;
info_sabcond.libprefix = libprefix;
info_sabcond.libprefix_Alib = libprefix_Alib;
info_sabcond.optLibs = optLibs;
info_sabcond.opticelib = opticelib;
info_sabcond.optInterpid = optInterpid;
info_sabcond.bands_opt = bands_opt;
info_sabcond.cntRmvl = cntRmvl;
info_sabcond.cntRmvl_ice = cntRmvl_ice;
info_sabcond.basenameWA = TRRIFdata.basenamesCDR.WA;
info_sabcond.wa = wa;
info_sabcond.bands = bands;
% info_sabcond.lines = lines;
info_sabcond.T_est = T_est;
info_sabcond.basename = basename_cr;
info_sabcond.save_dir = save_dir;


end

