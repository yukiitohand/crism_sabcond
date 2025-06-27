function [Asurficelib,infoAsurficelib] = ...
    crmsab_load_surficelibconv_wrapper(opt,wabasename,c)

switch opt
    case 1
        [A,infoA] = crmsab_load_lib_base( ...
            'Hapke_H2Oicelib_Grundy1998', 1, wabasename, c, ...
            'METHOD', 'interpCRISMspc', 'retainRatio',0.1);
        Asurficelib = A;
        infoAsurficelib = infoA;
    
    case 6
        [A,infoA] = crmsab_load_lib_base( ...
            'Hapke_H2Oicelib_Grundy1998', 6, wabasename, c, ...
            'METHOD','interpCRISMspc','retainRatio',0.1);
        Asurficelib = A;
        infoAsurficelib = infoA;
    case 7
        [A,infoA] = crmsab_load_lib_base( ...
            'Hapke_H2Oicelib_Grundy1998', 7, wabasename, c, ...
            'METHOD','interpCRISMspc','retainRatio',0.1);
        Asurficelib = A;
        infoAsurficelib = infoA;
    case 8
        [A,infoA] = crmsab_load_lib_base( ...
            'Hapke_H2Oicelib_Grundy1998', 8, wabasename, c, ...
            'METHOD','interpCRISMspc','retainRatio',0.1);
        Asurficelib = A;
        infoAsurficelib = infoA;
    case 9
        [A,infoA] = crmsab_load_lib_base( ...
            'Hapke_H2Oicelib_Grundy1998', 9, wabasename, c, ...
            'METHOD','interpCRISMspc','retainRatio',0.1);
        Asurficelib = A;
        infoAsurficelib = infoA;
    otherwise
        error('Opt %d is not defined',opt);
end





end