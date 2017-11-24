% updata beta
function[beta]=updatebeta(index,h,WW,Pdist,data,X,Init,A)
    if size(Pdist,1)>1
        Dis=Pdist(:,index);
    else
        Dis=GetDist(Pdist,index);
    end
    w=WW(Dis<=h,index);
    res=cell2mat(Init(Dis<=h,2));
    Y=data(:,Dis<=h)*(w./res);
    beta=(A{index,:}\X'*Y)';     
end