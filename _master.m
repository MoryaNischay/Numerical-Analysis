


=====BISECTION=====
a=input('Enter the value of a');
b=input('Enter the value of b');
tol=input('Enter the value of tolerance');
f=@(x) x^2-29;
if(f(a)*f(b)>0)
    disp('roots doess not lie in interval');
else
    mid=(a+b)/2;
    err=abs(a-b);
    while(err>tol)
        if((f(a)*f(mid))<0)
            b=mid;
        else
            a=mid;
        end
        mid=(a+b)/2;
        err=abs(a-b);
    end
end
disp(mid);




======FIXEDPOINT======
clc;
clear all;
f=@(x) atan(4*x);
a=input('Enter the starting value');
n=input('Enter the no of iterations');
tol=input('Enter the value of tolerance');
i=1;
x=f(a);
if(x==a)
    disp('root is a');
else
    while(i<=n)
        x=f(a);
        if(abs(x-a)<=tol)
            break;
        end
        a=x;
        i=i+1;
    end
end
disp(x);





====NEWTON=====
syms x;
z=@(x) cos(x)-x*exp(x);  
x0=input('Enter the value of a');
n=input('Enter the no of iterations');
tol=input('Enter the value of tolerance');
dif=diff(sym(z));
d=inline(dif);
i=1;
while(i<=n)
    x1=x0-((z(x0)/d(x0)));
    err=abs((x1-x0));
    if(err<=tol)
        disp(x1);
        break
    end
    x0=x1;
    i=i+1;
end





=====SECANT=====
f=@(x) x^2-17;
x0=input('Enter the value of a');
x1=input('Enter the value of b');
n=input('Enter the no of iterations');
tol=input('Enter the value of tolerance');
i=1;
while(i<=n)
    x2=x1-((x1-x0)/(f(x1)-f(x0)))*f(x1);
    err=abs(x2-x1);
    if(err<=tol)
        disp(x2);
        break
    end
    i=i+1;
    x0=x1;
    x1=x2;
end




=====GAUSSELE====
clc
A=[10 8 -3 1 16; 2 10 1 -4 9; 3 -4 10 1 10; 2 2 -3 10 11];
n=size(A,1);
for i=1:n-1
    for j=i+1:n
        A(j,:)=A(j,:)-(A(j,i)/A(i,i))*A(i,:);
    end
end
x(n)=A(n,n+1)/A(n,n);
for i=n-1:-1:1
    sum=0;
    for j=(i+1):n
        sum=sum+A(i,j)*x(j);
    end
    
    x(i)=(A(i,n+1)-sum)/A(i,i);
end
display(x);




=====LU======
clc
A = [2,-1,1;3,3,9;3,3,5];
n = size(A, 1);
L = eye(n);
for k = 1 : n
    L(k + 1 : n, k) = A(k + 1 : n, k) / A(k, k);
    for l = k + 1 : n
        A(l, :) = A(l, :) - L(l, k) * A(k, :);
    end
end
U = A;
U
L




====GAUSSSEIDEL=====
a=input('enter matrix a');
b=input('enter matrix b');
p=input('enter matrix c');
N=input('enter no of iterations');
n=length(b);
x=zeros(n,1);
y=p;
for j=1:N
    for i=1:n
        x(i)=(b(i)-a(i,[1:i-1,i+1:n])*p([1:i-1,i+1:n]))/a(i,i);
        p(i)=x(i);  
    end
      
        if(abs(x-y)<=0.0001)
            break
        end
        y=x;
end
disp(x);




======POWER=====
a=[4,1,0;1 20 1;0 1 4];
b=[1,1,1];
c = transpose(b)
n=input('enter no of iterations');
tol=input('Enter the value of tolerance');
for i=1:n
    v=a*c;
    x=max(abs(v));
    v1=v/x;
    if(abs(v1-c)<=tol)
        break;
    end
    c=v1;
end
disp(v1);
disp(x);



====LAGrange====
clc
clear
x=[1950 1960 1970 1980 1990 2000]
f=[151326 179323 203302 226542 249633 281422]
p=1965;
n=6;
for i=1:n
    l(i)=1.0;
    for j=1:n
        if j~=i;
            l(i)=((p-x(j))/(x(i)-x(j)))*l(i);
            
        end
    end 
end 
    
 
sum=0.0;
for i=1:n
    sum=sum+(l(i)*f(i));
end 
sum



====NEWTONDIVIDEd=====

n=4;
x=[1 1.5 2 2.5];
y=[2.7183 4.4817 7.3891 12.1825];
p=2.25;
for i=1:n
    D(i,1)=y(i);
end
for j=2:n
    for i=j:n
        D(i,j)=(D(i,j-1)-D(i-1,j-1))/(x(i)-x(i-j+1));
    end
end
disp(D);
prod = 1;
sum=D(1,1);
for i=2:n
    prod=prod*(p-x(i-1));
    sum=sum+prod*D(i,i);
end
disp(sum)





=====TRAPEZOIDAL====
f=@(x) (cos(x))^2;
a=input('Enter the value of a');
b=input('Enter the value of b');
n=input('Enter the no of subintervals');
h=(b-a)/n;
sum=0;
for i=1:n-1
    x=a+h*i;
    sum=sum+2*f(x);
end
sum=sum+f(a)+f(b);
sum=sum*(h/2);
disp(sum);




=====SIMSON====
f=@(x) (cos(x))^2;
a=input('Enter the value of a');
b=input('Enter the value of b');
n=input('Enter the no of subintervals');
h=(b-a)/n;
sum=0;
for i=1:n-1
    x(i)=a+h*i;
    y(i)=f(x(i));
    if mod(i,2)==0
        sum=sum+2*y(i);
    else
        sum=sum+4*y(i);
    end
end
sum=sum+f(a)+f(b);
sum=sum*(h/3);
disp(sum);



=====EULER====

f=@(x,y) -y + 2*cos(x);
a = 0;
b = 1;
h = 0.2;
n = (b - a) / h;
x = zeros(1,n);
y = zeros(1,n);
x(1) = a;
y(1) = 1;
for i = 1 : n
    x(i+1) = x(i) + h;
    y(i+1) = y(i) + h * f(x(i),y(i));
    y(i+1) = y(i) + (h/2) * (f(x(i+1), y(i+1)) + f(x(i),y(i)));
end
x
y


======RK======
clear all
f=@(t,y) -y+2*cos(t);
value=0;
final=1;
length=0.2;
y0=1;
while(value<final)
    k1=length*f(value,y0);
    k2=length*f(value+length/2,y0+k1/2);
    k3=length*f(value+length/2,y0+k2/2);
    value=value+length;
    k4=length*f(value+length,y0+k3);
    y0=y0+1/6*(k1+2*k2+2*k3+k4);
    disp(y0)
end

