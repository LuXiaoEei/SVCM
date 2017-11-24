function Dist=GetDist(Pdist,index)
    M=(1+sqrt(1+8*size(Pdist,2)))/2;
    Dist=zeros(M,1);
    for i=1:M
        if i~=index
            I=min(i,index);
            J=max(i,index);
            Dist(i,1)=Pdist((I-1)*(M-I/2)+J-I); 
        end
    end
end