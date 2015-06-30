function Lnorm = normalization_factors(Mnotes)
Lnorm = zeros(size(Mnotes, 1), size(Mnotes, 1));
for i=1:size(Mnotes,1)
    nnz = sum(Mnotes(i,:)~=0);
    if (nnz==0)
        Lnorm(i,i)=0;
    else
        Lnorm(i,i)=nnz^(-1);
    end
    
end
end