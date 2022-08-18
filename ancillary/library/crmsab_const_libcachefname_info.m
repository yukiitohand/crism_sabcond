function [infolibcache_fname] = crmsab_const_libcachefname_info(libname,opt)
infolibcache_fname = sprintf('%s%d_info.mat',libname,opt);
end