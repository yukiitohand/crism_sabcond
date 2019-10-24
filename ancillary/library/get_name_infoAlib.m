function [names] = get_name_infoAlib(infoA)
    names = cell(1,length(infoA));
    for i=1:length(infoA)
        switch infoA(i).spclib
            case 'CRISM spectral library'
                names{i} = infoA(i).spc_name;
            case 'RELAB spectral library'
                names{i} = [infoA(i).subType ' ' infoA(i).sampleName];
            case 'USGS splib06a'
                names{i} = infoA(i).name;
            case 'CRISM Type Spectral Library'
                names{i} = infoA(i).name;
            otherwise
                names{i} = '';
        end
    end
end