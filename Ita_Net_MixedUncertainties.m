%% Chen Automatic
clc
clear
N=50;
A=load('A.mat').A;
D=diag(sum(A,2));
L=D-A;

%% NE calculation
for i=1:N
    a_i=50;
    c=5;
    d=8;
    b_i=50-0.2*i;
    for j=1:N
        M(i,j)=c;
    end
    M(i,i)=2*a_i+2*c;
    q(i,1)=2*a_i*b_i-d;
end
NE=inv(M)*q;


X=10*ones(N,N);
sigma=5;
T_end=10000;
for k=1:T_end
    c_k=0.1/(log(k+1)+1+0.02*k^(3/4));
    c(k)=c_k;
    %c_k=0.05;
    %alpha_k=0.05/(log(k+1)+1+0.05*k^(3/4));
    alpha_k=c_k^(1.2);
    alpha(k)=alpha_k;
    %alpha_k=0.005;
    Omega=normrnd(0,0.1,1,N);
    %Omega=zeros(1,N);
    Y=X+kron(ones(N,1),Omega);
    for i=1:N
        x=diag(X);
        x_i=x(i);
        consen=0;
        PI=normrnd(0,1,N,N);
        for j=1:N
            consen=consen+A(i,j)*(Y(i,j)-Y(i,i))+A(i,j)*PI(i,j)*(Y(i,j)-Y(i,i));
        end
        x_i=x_i-alpha_k*gradient_f(X(:,i),i)+c_k*consen;
        %% estimate
        Y_others=Y;
        Y_others(i,:)=[];
        xi_others=Y_others(:,i);
        iSum_others=0;
        for j=1:N
            iSum_others=iSum_others+A(i,j)*(1+PI(i,j))*(Y_others(:,j)-Y_others(:,i));
        end
        xi_others=xi_others+c_k*iSum_others;

        if i==1
            X_i=[x_i;xi_others];
        elseif i==N
            X_i=[xi_others;x_i];
        else
            X_i(1:i-1)=xi_others(1:i-1);
            X_i(i)=x_i;
            X_i(i+1:end)=xi_others(i:end);
        end

        %存储起来
        X_STemp(:,i)=X_i;
    end
    X=X_STemp;
    %    x_star(:,k)=NE;
    S_x(:,k)=diag(X_STemp);
    Error_x(:,k)=abs(S_x(:,k)-NE);
    S_X(N*k-N+1:N*k,:)=X_STemp;
    T(k)=k;
end
x0=S_x;

%% Figure

figure (1)
for i=1:N
    plot(T,x0(i,:),'linewidth',1.2);
    hold on
end
axis([0,inf,-10,50]);

figure (2)
for i=1:N
    plot(T,Error_x(i,:),'linewidth',1.2);
    hold on
end
axis([0,inf,0,50]);

figure (3)
p1=plot(c);
hold on
p2=plot(alpha);
h1=legend([p1,p2],'C','Alpha');
set(h1,'fontsize',10)