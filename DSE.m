%double-stage estimation;
function [theta0,theta1,sigma,Vmu,V]=DSE(T1,C1,T2,C2,T3,C3,x1,x2,x3,stheta1,c)

T2t=T3*exp(stheta1*(x2-x3));
T1t=T3*exp(stheta1*(x1-x3));

t2=[T2,T2t];
t1=[T1,T1t];
t3=T3;

c3=C3;
c1=[C1,C3];
c2=[C2,C3];

x=[1,x1;1,x2;1,x3];

para1=mle(t1,'distribution','Weibull','censoring',c1);
Cov1=acov(log(para1(1,1)),1/para1(1,2),t1,c1);
para2=mle(t2,'distribution','Weibull','censoring',c2);
Cov2=acov(log(para2(1,1)),1/para2(1,2),t2,c2);
para3=mle(t3,'distribution','Weibull','censoring',c3);
Cov3=acov(log(para3(1,1)),1/para3(1,2),t3,c3);

y=[log(para1(1,1));log(para2(1,1));log(para3(1,1))];
Vmu=[Cov1(1,1),0,0;0,Cov2(1,1),0;0,0,Cov3(1,1)];
p=inv(x'*inv(Vmu)*x)*(x'*inv(Vmu)*y);
theta0=p(1);
theta1=p(2);
sigma=1/( 1/Cov1(2,2) +1/Cov2(2,2) +1/Cov3(2,2) ) * ( para1(1,2)^(-1)/Cov1(2,2) +para2(1,2)^(-1)/Cov2(2,2) +para3(1,2)^(-1)/Cov3(2,2) );

V=inv(x'*inv(Vmu)*x);
V(1,3)=0;
V(2,3)=0;
V(3,1)=0;
V(3,2)=0;
V(3,3)=1/( 1/Cov1(2,2) +1/Cov2(2,2) +1/Cov3(2,2) );

