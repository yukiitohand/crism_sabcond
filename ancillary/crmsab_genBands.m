function [ bands ] = crmsab_genBands( bands_opt )
% [ bands ] = crmsab_genBands( bands_opt )
%   generate "bands" used for processing from its option ID.

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


end

