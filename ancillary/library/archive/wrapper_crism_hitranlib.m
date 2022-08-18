function [Ahitranlib,infoAhitranlib] = wrapper_crism_hitranlib(opt,wabasename,c)


[Ah2o_hitran,infoh2o_hitran] = crmsab_load_lib_base('absxsecH2Olib_HITRAN',3,wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
[Aco2_hitran,infoco2_hitran] = crmsab_load_lib_base('absxsecCO2lib_HITRAN',3,wabasename,...
                                  c,'METHOD','interpCRISMspc','retainRatio',0.1);
switch opt
    case 1
        Ahitranlib = -Ah2o_hitran(:,7:34)*10^19*0.5;
        infoAhitranlib = infoh2o_hitran(7:34); 
    case 2
        Ahitranlib = cat(2,-Aco2_hitran(:,7:34),-Ah2o_hitran(:,7:34)*10^19*0.5);
        infoAhitranlib = merge_struct(infoco2_hitran(7:34),infoh2o_hitran(7:34)); 
    otherwise
        error('Opt %d is not defined',opt);
end

end