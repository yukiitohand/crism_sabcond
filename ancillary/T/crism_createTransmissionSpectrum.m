function [T_interpd] = crism_createTransmissionSpectrum(TRRIFdata,opt_sabcond,varargin)

% bands_opt = 4;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
%             case 'BANDS_OPT'
%                 bands_opt = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

% bands = crmsab_genBands(bands_opt);

% crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L');
% crism_obs.load_data(crism_obs.info.basenameIF,crism_obs.info.dir_trdr,'if');
% crism_obs.load_default();
% TRRIFdata = crism_obs.data.if;
TRRIFdata.load_basenamesCDR();
TRRIFdata.readWAi();

wa = TRRIFdata.wa;

data_sabcond = DATA_sabcond(TRRIFdata,opt_sabcond{:});

% load estimated T
load(data_sabcond.info.fname_supple, 'T_est');

[T_interpd] = crism_interp_TransSpc(T_est,data_sabcond.info.bands,TRRIFdata);

end


