%update cov
function[cov]=updatecov(index,espn,A,X)
    temp=A{index,1};
    cov=temp\X'*diag(espn(index,:).^2)*X/temp;       
end