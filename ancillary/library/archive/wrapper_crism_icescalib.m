function [Aicescalib,infoAicescalib] = wrapper_crism_icescalib(opt,wabasename,c)

switch opt
    case 4
        [AGrundy,infoAGrundy] = crmsab_load_lib_base('Q_H2Oicelib_Grundy1998',2,wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
        Aicescalib = AGrundy;
        infoAicescalib = infoAGrundy;
    otherwise
        error('Opt %d is not defined',opt);
end





end