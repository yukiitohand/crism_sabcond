function [optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib] = decompose_libprefix(libprefix)

switch length(libprefix)
    case 7
        optCRISMspclib  = str2num(libprefix(4));
        optRELAB        = str2num(libprefix(5));
        optUSGSsplib    = str2num(libprefix(6));
        optCRISMTypeLib = str2num(libprefix(7));
        opticelib = '';
    case 8
        optCRISMspclib  = str2num(libprefix(4));
        optRELAB        = str2num(libprefix(5));
        optUSGSsplib    = str2num(libprefix(6));
        optCRISMTypeLib = str2num(libprefix(7));
        opticelib       = str2num(libprefix(8));
    otherwise
        error('libprefix %s is invalid',libprefix);
end