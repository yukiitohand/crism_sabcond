function [wa_identfr] = crmsab_get_libWAIdentifier(wabasename)
% [wa_identfr] = crmsab_get_libWAIdentifier(wabasename)
% Get identfier for CDR WA of the given basename
% basically sclk is set to zerom since files only different in sclk are
% identical.
% 'CDR410803692813_WA0000000L_3' --> 'CDR410000000000_WA0000000L_3'


propWA = crism_getProp_basenameCDR4(wabasename);
propWA.sclk = 0;

wa_identfr = crism_get_basenameCDR4_fromProp(propWA);

end