% 计算加权最小二乘的权重 核函数:(1-u)+
function Wh = GetWh(Pdist_xyz,h,index)
    Pdist_x=Pdist_xyz(:,1);
    Pdist_y=Pdist_xyz(:,2);
    Pdist_z=Pdist_xyz(:,3);
    u1=GetDist(Pdist_x,index)/h;
    u2=GetDist(Pdist_y,index)/h;
    u3=GetDist(Pdist_z,index)/h;
%     u1=GetDist(pdist(postion(:,1)),index)/h;
%     u2=GetDist(pdist(postion(:,2)),index)/h;
%     u3=GetDist(pdist(postion(:,3)),index)/h;
    Wh=(1-u1).*(1-u1>0).*(1-u2).*(1-u2>0).*(1-u3).*(1-u3>0)/h^3;
end

    
