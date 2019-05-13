function [ trans_spcs ] = load_T_sclk_closest( sclk,varargin )
% [ trans_spcs ] = load_T_sclk_closest( sclk,varargin )
% transmission is loaded based on the closest in terms of the input sclk.
% The transmission files to be read
% Optional parameters
%  'ADDITIONAL_SUFFIX': ad
%      (default) 'lam05'


additional_suffix = 'lam05';
% binning = ''; wv_filter = '';
% overwrite = 0;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'ADDITIONAL_SUFFIX'
                additional_suffix = varargin{i+1};
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

%% select the closest FFC

[ dirname_ffc,dirname_ffc_mcb,dirname_ffc_mca ] = select_closestFFC(sclk);


%% resolve the filenames
T_matfname = sprintf('TList_%s_%s.mat',dirname_ffc,additional_suffix);
if ~exist(T_matfname,'file')
    error('%s does not exist.',T_matfname);
end

load(T_matfname,'TList');
trans_spcs = TList;

end