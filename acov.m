function [COV]=acov(m,s,t,C)

n=length(t);
y=log(t);


%%%%%%%%%%%%%%%%%% derivative with mu and mu %%%%%%%%%%%%%%%%%%%%%%
l=0;
for i=1:1:n
    %if t1(i)<c
        l=l-exp( (y(i)-m)/s ) / s^2;
    %else
    %    l=l-exp( (y(i)-m)/s ) / s^2;
    %end
end

f11=l;

%%%%%%%%%%%%%%%%%% derivative with mu and sigma%%%%%%%%%%%%%%%%%%%%%%
l=0;
for i=1:1:n
    if C(i)==0
        l=l + 1/s^2 - exp( (y(i)-m)/s ) / s^2 - exp( (y(i)-m)/s )*(y(i)-m)/s /s^2;
    else
        l=l- exp( (y(i)-m)/s ) / s^2 - exp( (y(i)-m)/s )*(y(i)-m)/s /s^2;
    end
end

f12=l;

%%%%%%%%%%%%%%%%%% derivative with sigma and sigma%%%%%%%%%%%%%%%%%%%%%%
l=0;
for i=1:1:n
    if C(i)==0
        l=l + 1/s^2 + 2*(y(i)-m)/s/s^2 - 2*exp( (y(i)-m)/s )*(y(i)-m)/s / s^2 - exp( (y(i)-m)/s )*((y(i)-m)/s)^2 /s^2;
    else
        l=l - 2*exp( (y(i)-m)/s )*(y(i)-m)/s / s^2 - exp( (y(i)-m)/s )*((y(i)-m)/s)^2 /s^2;
    end
end

f22=l;

F2=-[f11,f12;f12,f22];
COV=inv(F2);