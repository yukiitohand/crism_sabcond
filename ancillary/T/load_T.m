function [ trans_spcs ] = load_T( varargin )
% [ trans_spcs ] = load_ADR_VS( varargin )
%
% Optional parameters
%    'BINNING': binning mode {0,1,2,3,''}
%               (default) ''
%    'WAVELENGH_FILTER': wavelenth filter {0,1,2,3,''}
%                       (default) ''
%    'VERSION: (default) 8
%    'OBS_ID_SHORT': (defaut) ''
%    'T_MODE' : mode for transmission selection,
%                1: single manual selection, need to specify "ADR_BASE".
%                2: 
%                3: concatenatation of all the available transmission data
%               (default) 2
%    'ADR_BASE': basename of the manually selected data, used only with
%                T_MODE = 1.  (default) []
%    'ARTIFACT': specify how to deal with artifact. {'subtraction', 'none'}
%                default 'subtraction'
%    'ARTIFACT_IDX': specify which data is used as artifact {2,3}. 2: new
%                    Patrick C. McGuire method and 3: old Pelky method.
%                    (default) 2
%    'OVERWRITE' : whether or not to overwrite the cache file
%                  (default) 0
%
global localCRISM_PDSrootDir

t_mode = 2;
adr_base = '';
artifact_idx = 2;
opt_artifact = 'subtraction';
binning = ''; wv_filter = ''; vr = '8'; obs_id_short = '';
overwrite = 0;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING'
                binning = varargin{i+1};
                if isnumeric(binning)
                    binning = sprintf('%1d',binning);
                end
            case 'WAVELENGTH_FILTER'
                wv_filter = varargin{i+1};
                 if isnumeric(wv_filter)
                    wv_filter = sprintf('%1d',wv_filter);
                 end
            case 'VERSION'
                vr = varargin{i+1};
                if isnumeric(vr)
                    vr = sprintf('%1d',vr);
                end
            case 'OBS_ID_SHORT'
                obs_id_short = varargin{i+1};
            case 'T_MODE'
                t_mode = varargin{i+1};
            case 'ADR_BASE'
                adr_base = varargin{i+1};
            case 'ARTIFACT'
                opt_artifact = varargin{i+1};
            case 'ARTIFACT_IDX'
                artifact_idx = varargin{i+1};
            case 'OVERWRITE'
                overwrite = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

if isempty(adr_base)
    propADRVSPtr = create_propADRVSbasename();
    if ~isempty(binning), propADRVSPtr.binning = binning; end
    if ~isempty(wv_filter), propADRVSPtr.wavelength_filter = wv_filter; end
    if ~isempty(vr), propADRVSPtr.version = vr; end
    if ~isempty(obs_id_short), propADRVSPtr.obs_id_short = obs_id_short; end
else
    propADRVSPtr = getProp_basenameADRVS(adr_base);
end

%%
load('TList.mat','TList');
trans_spcs = TList;

end