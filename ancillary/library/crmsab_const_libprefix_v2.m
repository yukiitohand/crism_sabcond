function libprefix = crmsab_const_libprefix_v2(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib,opthitranlib)

libprefix = sprintf('Lib%d%d%d%d%d%d',optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib,opthitranlib);

end