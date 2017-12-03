%ͨ��WARM�������ƾ���
%input data��ͼ�����ݼ�n*p ch�����봰����� Cn�����ƴ��� Kloc���Ƿ�ʹ�þ���� Kst���Ƿ�ʹ�����ƺ�S0��ֹͣ��⿪ʼ�������� S�������ܴ��� Dis��������� X���Ա���
%output ���ƾ���WW ÿi�б�ʾ��i���������������ص����ƶ���
function[Init0]=WARM(data,ch,Cn,S0,S,Pdist,X,Init)
	Dis=squareform(Pdist);
    %Pdist=squareform(Pdist);% ###A
    fprintf('----------Initializing-----%s----------\n',datestr(now()));
    fprintf('----------Initializing Beta-----%s----------\n',datestr(now()));
    Init0=Init;
    fprintf('----------Iteration Begin -----%s----------\n',datestr(now()));
    dict1=false(size(data,2),1);%�������ֹͣ�ı��
    for i=1:S
        fprintf('----------��ʼ�ص�%d�ε��� ----------\n',i);
        fprintf('------------Update WW -----%s----------\n',datestr(now()));      
        try
            W_loc=((1-Dis/ch^i).*((1-Dis/ch^i)>=0)); 
        catch
            W_loc=(1-Pdist/ch^i).*((1-Pdist/ch^i)>=0);
        end
        fprintf('------------Update WW Matrix -----%s----------\n',datestr(now()));
        tic;
        WW=Weight(Init,W_loc,Cn);
        toc;
%         
%         disp(sum(sum(WW-WWW)));
        
        fprintf('------------Update A -----%s----------\n',datestr(now()));
        A=cell(size(data,2),1);
		tic;
        for index=1:size(data,2)
            A{index,:}=CompA(index,ch^i,WW,Pdist,Init,X);
        end
		toc;
		
        fprintf('------------Update Beta -----%s----------\n',datestr(now()));
        beta=cell(size(data,2),1);

		tic;
        for index=1:size(data,2)
           beta{index,:}=updatebeta(index,ch^i,WW,Pdist,data,X,Init,A);
        end
		toc;
		
        beta(dict1,:)=Init0(dict1,1);%ֹͣ�Ĳ����е���
        if i>S0
%                ֹͣ׼�� 
            D=stopctiter(Init,Init3);
            dict0=D>chi2inv(0.8,2);
            dict=(~dict1)&dict0;%�������Ҫֹͣ�ı��
            dict1=dict1|dict0;%�������ֹͣ�ı��
            beta(dict,:)=Init0(dict,1);%ֹͣ                  
        end
        Init(:,1)=beta;
        fprintf('------------Update Cov -----%s----------\n',datestr(now()));
        espn=zeros(size(data,2),size(data,1));
		tic;
        for index=1:size(data,2)
            espn(index,:)=epsn(index,WW,Pdist,ch^i,Init,data,X);
        end
		toc;
		
        tic;
        parfor index=1:size(data,2)
            Init{index,3}=updatecov(index,espn,A,X);
        end
		toc;	
        if i==S0
            Init3=Init; 
        end
        Init0=Init;
    end
    fprintf('----------Iteration End -----%s----------\n',datestr(now())); 
end