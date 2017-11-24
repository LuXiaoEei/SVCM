cd E:\mypackages\SVCM
clear;
postion = [repmat(reshape(repmat(1:64,64,1),1,64*64),1,8);
    repmat(1:64,1,64*8);
    reshape(repmat(1:8,64*64,1),1,64*64*8)]';
postion = [repmat(reshape(repmat(1:20,20,1),1,400),1,10);
    repmat(1:20,1,200);
    reshape(repmat(1:10,400,1),1,4000)]';
% % Dis=squareform(Pdist);
%% 生成数据
[Beta,A,Fai,Eta,X,Y_norm,Y_chip]=GenerateData(postion,80);
Pdist=pdist(postion);
% GetDist(Pdist,1);
%% 选择数据集,进行SVCM
Ydata=Y_norm;
%% stage1
% estimate beta0
beta0=zeros(size(Ydata,2),size(X,2));
% 无截距项
%  X=[ones(length(X),1),X];
oumu=inv(X'*X);
for index=1:size(Ydata,2)
   beta0(index,:)=(oumu*X'*Ydata(:,index))'; 
end
% estimate eta
res=Ydata-X*beta0'; % residual
Res=sum(res.^2,1)'/(size(Ydata,1)-size(X,2));
Init=cell(size(Ydata,2),3);  
for i=1:size(Init,1)
    Init(i,:)={beta0(i,:) Res(i,:) oumu*Res(i,:)};
end
clear oumu Res;
% h=10*size(Ydata,2)^-(1/5); %windows width
h=5;

tic;
Eta_est=ones(size(Ydata));
for index=1:size(Ydata,1)
%     Pdist_x=pdist(postion(:,1));
%     Pdist_y=pdist(postion(:,2));
%     Pdist_z=pdist(postion(:,3));
%     Pdist_xyz=[Pdist_x,Pdist_y,Pdist_z];
%     Wh = GetWh(Pdist_xyz,h,ii)
    temp=zeros(size(Ydata,2),1);
     tic;
    for ii=1:size(Ydata,2)
%         u=Dis(:,ii); % ###A
        u=GetDist(Pdist,ii);
        Wh=(1-u/h).*(u<h);
        Zh=[ones(size(Ydata,2),1),postion(:,1)-postion(ii,1),postion(:,2)-postion(ii,2),postion(:,3)-postion(ii,3)]; 
%         temp(ii,1)=[1,0,0,0]*((Zh'*diag(Wh)*Zh)\Zh'*diag(Wh)*res(index,:)');       
        temp(ii,1)=[1,0,0,0]/(Zh'*(repmat(Wh,1,4).*Zh))*Zh'*(Wh.*res(index,:)');%Zh不是4000*4（3维）则需要更改代码             
    end
    Eta_est(index,:)=temp;
    toc;
end
toc;

%estimate Fai
% Eta_norm=Eta_est-ones(size(Eta_est,1),1)*mean(Eta_est);
K=svd(Eta_est);
[U,~,V]=svds(Eta_est,5);

% K=svd(cov_eta);
% [U,~,V]=svds(cov_eta,5);

num=1;
siz=20;
figure(1)
imagesc(reshape(V(postion(:,3)==6,num),siz,siz))
figure(2)
v=Fai{1,:}(:,num);
imagesc(reshape(v(postion(:,3)==1,1),siz,siz))

save D:\共享文件\result_4000h5.mat

% 样本协方差
res2=res-Eta_est;
cov_esp=mean(res2.^2,1)';% 随机误差的方差
cov_eta=Eta_est'*Eta_est/(size(X,1)-size(X,2));%随机效应的协方差阵,当内存不够的时候，可以使用CovEta计算第index列的协方差.

%% stage2
% MASS
ch=1.1;
Cn = size(Ydata,1)^0.4 * chi2inv(0.8,1);
S0=3;
S=10;
warm=WARM(Ydata,ch,Cn,S0,S,Pdist,X,Init);
%estimate
beta_final=cell2mat(warm(:,1));
% imagesc(reshape(beta_final(postion(:,3)==1,2),20,20))
% imagesc(reshape(beta0(postion(:,3)==1,2),20,20))
% imagesc(reshape(Beta(postion(:,3)==1,2),20,20))
