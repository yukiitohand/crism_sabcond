function [CRISMspclib,libs_CRISMspclib] = load_CRISMspclib(opt_CRISMspclib)
    % crism spectral library
    switch opt_CRISMspclib
        case 1
            % mineral + rocks
            [CRISMspclib] = readCRISMspclib();
            libs_CRISMspclib = {'carbonate','inosil','nesosil','nitrates',...
                'oxide','phosphate','phylosil','sorosil','sulfate',...
                'tectosil','rocks'};
        case 2
            % all
            [CRISMspclib] = readCRISMspclib();
            libs_CRISMspclib = 'all';
        case 3
            % mineral
            [CRISMspclib] = readCRISMspclib();
            libs_CRISMspclib = {'carbonate','inosil','nesosil','nitrates',...
                'oxide','phosphate','phylosil','sorosil','sulfate',...
                'tectosil'};
        otherwise
            error('opt_CRISMspclib %d is not defined',opt_CRISMspclib);
    end
end