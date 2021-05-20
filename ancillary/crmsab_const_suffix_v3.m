function [suffix] = const_suffix_v3(mt,additional_suffix)


if isempty(additional_suffix)
    suffix = ['_' mt];
else
    suffix = ['_' mt '_' additional_suffix];
end

end