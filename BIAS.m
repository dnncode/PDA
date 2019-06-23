%double-stage estimation;
function [output]=BIAS(r1,r2,r3,x1,x2,x3,x0,e,sigma,V)

error1=exp(e*(x1-x3)/sigma);
error2=exp(e*(x2-x3)/sigma);
S=[1,x1;1,x2;1,x3];

% computation of bias
for i=1:1:1000
    X1=chi2rnd(2*(r1+r3));
    X2=chi2rnd(2*(r2+r3));
    X3=chi2rnd(2*r3);
    
    bias1(i)=log( X1/2 + (error1-1)*X3/2 ) - log(r1 + r3);
    bias2(i)=log( X2/2 + (error2-1)*X3/2 ) - log(r2 + r3);
    bias3(i)=log( X3/2 ) - log( r3);
end
Bias1=sigma*mean(bias1);
Bias2=sigma*mean(bias2);
Bias3=sigma*mean(bias3);

Bias=[Bias1;Bias2;Bias3];

output=[1,x0]*inv(S'*inv(V)*S)*(S'*inv(V)*Bias);


