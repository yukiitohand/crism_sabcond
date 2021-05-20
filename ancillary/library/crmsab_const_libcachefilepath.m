function [libcachefilepath] = crmsab_const_libcachefilepath(pdir,masterbase,c)
libcachefilepath = joinPath(pdir,sprintf('%s_c%03d.mat',masterbase,c));
end