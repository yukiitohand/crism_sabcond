function [libcachefname] = crmsab_const_libcachefname(libname,opt,wabasename,method,retainRatio,c)
switch method
    case 'interp1'
        libcachefname = sprintf('%s%d_%s_%s_c%03d',libname,opt,wabasename,method,c);
    case 'interpCRISMspc'
        libcachefname = sprintf('%s%d_%s_r%d_c%03d',libname,opt,wabasename,retainRatio*10,c);
    otherwise
        error('method %s is not valid',method);
end
end