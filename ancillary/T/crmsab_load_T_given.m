function [ trans_spcs ] = crmsab_load_T_given( obs_id,varargin )
% [ trans_spcs ] = crmsab_load_T_given( obs_id,varargin )
% transmission is loaded based on the give obs_id
% INPUTS
%  obs_id: string, 8characters or less
% Optional parameters
%  'ADDITIONAL_SUFFIX': ad
%      (default) 'lam05'
%  'DIRPATH': 


additional_suffix = 'lam05';
% binning = ''; wv_filter = '';
% overwrite = 0;
dirpath = '';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'ADDITIONAL_SUFFIX'
                additional_suffix = varargin{i+1};
            case 'DIRPATH'
                dirpath = varargin{i+1};
%             case 'BINNING'
%                 binning = varargin{i+1};
%                 if isnumeric(binning)
%                     binning = sprintf('%1d',binning);
%                 end   
%             
%             case 'WAVELENGTH_FILTER'
%                 wv_filter = varargin{i+1};
%                  if isnumeric(wv_filter)
%                     wv_filter = sprintf('%1d',wv_filter);
%                  end
%             case 'OVERWRITE'
%                 overwrite = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

%% get ffc basename
obs_id8c = sprintf('%08s',obs_id);


%% resolve the filenames
T_matfname = joinPath(dirpath,sprintf('TList_FFC%s_%s.mat',obs_id8c,additional_suffix));
if ~exist(T_matfname,'file')
    error('%s does not exist.',T_matfname);
end

load(T_matfname,'TList');
trans_spcs = TList;

end