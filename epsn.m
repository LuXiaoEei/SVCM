% ¼ÆËã¼ÓÈ¨²Ð²î
function[esp]=epsn(index,WW,Pdist,h,Init,data,X)
    if size(Pdist,1)>1
        Dis=Pdist(:,index);
    else
        Dis=GetDist(Pdist,index);
    end
    Beta=cell2mat(Init(Dis<=h,1))';
    Y=data(:,Dis<=h);
    w=WW(Dis<=h,index);
    res=cell2mat(Init(Dis<=h,2));
    esp=(w./res)'*(Y-X*Beta)';  %1by100  
end