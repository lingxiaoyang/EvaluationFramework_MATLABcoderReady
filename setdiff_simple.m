function c = setdiff_simple(a,b)
% ONLY SUPPORT TWO IN ARGS and ONE OUT ARG only FOR ROW VEC
%SETDIFF Set difference.
%   C = SETDIFF(A,B) for vectors A and B, returns the values in A that 
%   are not in B with no repetitions. C will be sorted.

    % Call ISMEMBER to determine list of non-matching elements of A.
    logUA = ~(ismember(a,b));
    c = a(logUA);
    
    % Call UNIQUE to remove duplicates from list of non-matches.
    c = unique(c);
end
