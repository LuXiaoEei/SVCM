% compute A
function[AA]=CompA(index,h,WW,Pdist,Init,X)
    if size(Pdist,1)>1
        Dis=Pdist(:,index);
    else
        Dis=GetDist(Pdist,index);
    end
    w=WW(Dis<=h,index);
    res=Init(Dis<=h,2);
    AA=(X'*X)*sum(w./cell2mat(res));
end