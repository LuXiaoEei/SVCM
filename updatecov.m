%update cov
function[cov]=updatecov(index,espn,A,X)
    temp=A{index,1};
    cov=temp\X'*(repmat((espn(index,:).^2)',1,size(X,2)).*X)/temp;       
end