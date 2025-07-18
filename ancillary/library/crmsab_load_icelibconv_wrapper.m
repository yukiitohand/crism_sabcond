function [Aicelib,infoAicelib] = crmsab_load_icelibconv_wrapper(opt,wabasename,c)

switch opt
    case 1
        [AGrundy,infoAGrundy] = crmsab_load_lib_base('abscoeffH2Oicelib_Grundy1998',[],wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
        Aicelib = -AGrundy*10;
        infoAicelib = infoAGrundy;
    case 2
        [AGrundy,infoAGrundy] = crmsab_load_lib_base('abscoeffH2Oicelib_Grundy1998',[],wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
        [AHansen,infoAHansen] = crmsab_load_lib_base('abscoeffCO2icelib_Hansen',[],wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
        Aicelib = [-AGrundy*10 -AHansen*10^2];
        infoAicelib = merge_struct(infoAGrundy,infoAHansen);
    case 3
        [AGrundy,infoAGrundy] = crmsab_load_lib_base('Q_H2Oicelib_Grundy1998',3,wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
        Aicelib = AGrundy;
        infoAicelib = infoAGrundy;
    case 4
        [AGrundy,infoAGrundy] = crmsab_load_lib_base('Q_H2Oicelib_Grundy1998',2,wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
        Aicelib = AGrundy;
        Aicelib = normalizevec(Aicelib,1);
        infoAicelib = infoAGrundy;
    case 5
        [AGrundy,infoAGrundy] = crmsab_load_lib_base('Q_H2Oicelib_Grundy1998',4,wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
        Aicelib = AGrundy;
        infoAicelib = infoAGrundy;
    case 6
        [AGrundy,infoAGrundy] = crmsab_load_lib_base('Q_H2Oicelib_Grundy1998',5,wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
        Aicelib = AGrundy;
        infoAicelib = infoAGrundy;

    case 9
        [AGrundy,infoAGrundy] = crmsab_load_lib_base('Q_H2Oicelib_Grundy1998',9,wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
        Aicelib = AGrundy;
        infoAicelib = infoAGrundy;

    otherwise
        error('Opt %d is not defined',opt);
end





end