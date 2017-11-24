% 计算第index点与其他点的协方差
function coveta =CovEta(Eta_est,p,index)
    n=size(Eta_est,1);
    coveta=(Eta_est(:,index)'*Eta_est/(n-p))'; 
end