function [Aall_cachefname] = crmsab_constAlibcachefname_Aall(libprefix,wabasename,optInterpid,c)

Aall_cachefname = sprintf('Alib%s_%s_interp%d_Aall_c%03d.mat',libprefix,wabasename,optInterpid,c);

end