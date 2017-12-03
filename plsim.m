function [g, theta, beta] = plsim(x, z, y, h)

% model: y = beta^T*z + g(theta^T*x) + e;
%   INPUT:  X: n x p;
%           Z: n x q; 
% 			   h: bandwidth 
%   OUTPUT: theta, beta, g.

[n,p]=size(x);
[n,q] = size(z);
pq = p+q;

onep = ones(p,1);
onen = ones(n,1);
B = eye(p);
nd = 1;

a = zeros(n,1);
h2 = 2*n^(-2/(p+4));
invzz = inv(z'*z)*z';
eyep1 = eye(p+1)/n^2;
% for iter = 1:p;
   ye = y - a;
   beta = invzz*ye;
   ye = ye - z*beta;
	ab = ones(p,n);
	for i = 1:n;
   	xi = x - repmat(x(i,:),n,1);
   	kernel = exp(-sum((xi*B).^2,2)/h2);
   	onexi = [onen xi];
   	xk = onexi.*repmat(kernel, 1, p+1);
      abi = inv(xk'*onexi+eyep1)*xk'*ye;
   	ab(:,i) = abi(2:p+1);
   	a(i) = abi(1);
	end;
% 	[B0 D] = eig(ab*ab');
% 	[D I] = sort(diag(D));
% 	B = B0;
% 	for k = 1:p;
%    	B0(:,k) = B(:,I(k));
% 	end
% 	for k = 1:nd;
%    	  B(:,k) = B0(:,p-k+1);
% 	end;  
%    B = B(:,1:max(1,p-iter));
% end;
% theta = B;
theta=mean(ab,2)/sqrt(sum(mean(ab,2).^2));

h2 = (2*h*h);
eyepq = eye(pq)/n^2;
for k = 1:20;
    tmp = repmat(x*theta, 1, n);
    d = (tmp-tmp').^2;
    ker = exp(-d/h2);
    ker = ker./repmat(sum(ker,2),1,n);
    md = zeros(n,p*p);
    mc = zeros(n,p);
    me = zeros(n,p);
    
    D22 = 0;
    C2 = 0;
    mcz = zeros(n,q);
    mez = zeros(n,q);
    for i = 1:n;
        tmp = x-repmat(x(i,:),n,1);
        tmp1 = repmat(ker(:,i),1,p).*tmp;
        md(i,:) = reshape(tmp1'*tmp,1,p*p);
        mc(i, :) = sum(tmp1,1);
        me(i, :) = y'*tmp1;
        
        z1 = repmat(ker(:,i),1,q).*z;
        D22 = D22 + z1'*z;
        mcz(i, :) = sum(z1,1);
        C2 = C2 + y'*z1;
        
        mdd(i,:) = reshape(tmp'*z1,1,q*p);
     end;
     
     for j = 1:5
        ye = y - z*beta;
        tmp = repmat(x*theta, 1, n);
        d1 = tmp-tmp';
        ker1 = d1.*ker;
        s2 = sum(d1.*ker1,2);
        s1 = sum(ker1, 2);
        s = sum(ker,2);
        d = s2.*s - s1.*s1 + 1/n^2;
        a = (ker*ye.*s2 - ker1*ye.*s1)./d;
        b = (ker1*ye.*s - ker*ye.*s1)./d;
        D = reshape((b.*b)'*md, p, p);
        C = b'*me - (b.*a)'*mc;
        
        Cz = C2 - a'*mcz; 
        D12 = reshape(b'*mdd, p, q);
        D = [D D12
           D12' D22];
        C = [C Cz];
        
        theta = inv(D+eyepq)*C';
        beta = theta(p+1:pq); theta = theta(1:p);
        theta = sign(theta(1))*theta/sqrt(theta'*theta);
    end
 end    
 
 g = a;
