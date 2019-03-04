function libprefix = const_libprefix(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib,opthitranlib)

libprefix = sprintf('Lib%d%d%d%d%d%d',optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib,opthitranlib);

end