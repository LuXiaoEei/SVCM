x=unifrnd(-1,1,100,6);
z=normrnd(0,1,100,5);
theta=[0.75 0.5 -0.25 -0.25 0.25 0.3];
theta=theta/sqrt(sum(theta.^2));
beta=[1 2 3 4 5];
g=60*exp(-(x*theta'-0.5).^2);
y=z*beta'+g+normrnd(0,0.6,100,1);
h=0.1;

[G,Theta,Beta] = plsim(x, z, y, h);
Y=z*Beta+G;

figure(1)
plot(-(x*theta'-0.5).^2,g,'.');
hold on;
plot(-(x*theta'-0.5).^2,G,'*')
hold off;

figure(2)
plot(Y,'DisplayName','Y');hold on;plot(y,'DisplayName','y');hold off;


