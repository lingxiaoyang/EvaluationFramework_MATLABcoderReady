function c = intersect_simple(a,b)
% ONLY SUPPORT TWO IN ARGS and ONE OUT ARG and ONLY COLUMN VEC
%INTERSECT Set intersection.
%   C = INTERSECT(A,B) for vectors A and B, returns the values common to
%   the two vectors with no repetitions. C will be sorted.

    % reverse a and b for better performance
    c = unique(b(ismember(b,a)));
    c = [c; a([])]; % make sure output has correct type
end
