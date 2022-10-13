function [ bands ] = crmsab_genBands_v2(wavelength_filter, bands_opt,binning,sclk )
% [ bands ] = crmsab_genBands_v2( wavelength_filter, bands_opt )
%   generate "bands" used for processing from its option ID.

switch wavelength_filter
    case 0
        switch bands_opt
            case 3
                bands = [4:97 105:244];
            case 4
                bands = [4:244]';
            case 5
                bands = [4:252];
            case 6
                bands = [4:250];
            case 7
                bands = [1:250];
            case 8
                bands = [1:252];
            otherwise
                error('Option %d is not defined yet');
        end
    case 2
        switch binning
            case 3
                switch sclk
                    case 947778566
                        switch bands_opt
                            case 4
                                bands = 1:135;
                            case 6
                                bands = 1:138;
                            otherwise
                                 error('Option %d is not defined yet');
                        end
                    otherwise
                        error('sclk=%d is not considered yet',sclk);
                end
            otherwise
                error('wavelength_filter=%d and binning=%d is not considered yet',wavelength_filter,binning);
        end
end


end