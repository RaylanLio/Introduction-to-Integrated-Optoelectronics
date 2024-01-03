clear all clc;

N=20;           %将光栅划分为N份
n=zeros(1,N+4); %n(1)->n0  
neff=linspace(3.291,3.449,1000);
L=length(neff); 
gamma=zeros(L,N+4); %gamma(1)->gamma0  
T=zeros(L,N+3);     %T(1)->T0
A=T;        %系数
B=T;
d=zeros(1,N+2);     %记录每一层的厚度
period=250e-9;      %周期
wavelength=880e-9;
k0=2*pi/wavelength; 
d(1)=0.1e-6;        %d(1)->d1
NEFF=[];
r=[];
b=[];

n(1)=3.29;      %0区折射率
n(2)=3.59;      %1区折射率
n(3)=3.45;      %2区折射率
n(end)=3.29;    %N+3区折射率
step_h=(0.1e-6-0.1e-6/1e6)/50;  %扫描的h的步长
A_temp=[];
B_temp=[];
h=0.1e-6/1e6:step_h:0.1e-6;     %扫描的h的范围0~0.1e-6
% % disp(['倒计时开始！']);
for p=1:length(h)
% %     disp(['还有',num2str(length(h)-p),'步完成！']);
    d(2)=0.15e-6-h(p);    %计算2区厚度
    h0=h(p)/N;  %光栅层均分，每层厚度一样
    d(3:end)=h0;       
    b=[b;d];        %将不同h下的每层厚度存在b中
    %计算平均折射率
    h_temp=h0*[1:20]-h0/2; %每一层的中间线距光栅顶的高度
    x_mid=period*h_temp/h(p);   %中间
    x2=period-x_mid;
    n(4:N+3)=sqrt((x_mid*3.29^2+x2*3.45^2)/period);%每一层的平均折射率
    err=[];
    for k=1:L
        gamma(k,:)=k0*sqrt(abs(neff(k)^2-n.^2));
        T(k,:)=gamma(k,1:end-1)./gamma(k,2:end);
        A(k,1)=1;         %计算A1
        B(k,1)=T(k,1);      %计算B1
        for m=2:N+2     %计算2~N+2区的系数
            if neff(k)<n(m)
                A(k,m)=A(k,m-1)*cos(gamma(k,m)*d(m-1))+...
                    B(k,m-1)*sin(gamma(k,m)*d(m-1));
                B(k,m)=T(k,m)*(B(k,m-1)*cos(gamma(k,m)*d(m-1))-...
                    A(k,m-1)*sin(gamma(k,m)*d(m-1))); 
            else
                A(k,m)=A(k,m-1)*cosh(gamma(k,m)*d(m-1))+...
                    B(k,m-1)*sinh(gamma(k,m)*d(m-1));
                B(k,m)=T(k,m)*(B(k,m-1)*cosh(gamma(k,m)*d(m-1))+...
                    A(k,m-1)*sinh(gamma(k,m)*d(m-1)));
            end
        end
        if neff(k)<n(N+3)   %计算N+3区的系数
            A(k,N+3)=A(k,N+2)*cos(gamma(k,N+3)*d(N+2))+...
                    B(k,N+2)*sin(gamma(k,N+3)*d(N+2));
        else
            A(k,N+3)=A(k,N+2)*cosh(gamma(k,N+3)*d(N+2))+...
                    B(k,N+2)*sinh(gamma(k,N+3)*d(N+2));
        end
        if neff(k)<n(N+2)   %计算不同neff下的本征方程值
        err_temp=tan(gamma(k,N+3)*h0)-(A(k,end-1)+...
            B(k,end-1)*T(k,end))/(A(k,end-1)*T(k,end)-B(k,end-1));
        else
            err_temp=tanh(gamma(k,N+3)*h0)+(A(k,end-1)+...
            B(k,end-1)*T(k,end))/(A(k,end-1)*T(k,end)+B(k,end-1));
        end
        err=[err err_temp];
    end
    err=abs(err);
    flag=find(err==min(err));   %寻找对应本征方程解的neff序数
    A_temp=[A_temp;A(flag,:)];  %将对应的场系数存起来
    B_temp=[B_temp;B(flag,:)];
    r=[r;gamma(flag,:)];        %将对应的gamma存起来
    NEFF=[NEFF neff(flag)];     
end
figure;
plot(h,NEFF),xlabel('h/m'),ylabel('n_e_f_f')
%title('n_e_f_f—h曲线');  

A=A_temp;
B=B_temp(:,1:end-1);
save('data.mat','A','B','r','b','h','n','NEFF','N');
clear all,close all,clc;
load data.mat;
period=250e-9;
beta0=2*pi/period;  
k0=2*pi/880e-9;
L=1001;              %每一层的取点数
Lh=length(h);       %选取的h的点数
%单独计算第0层的P
x=linspace(-1e-4,0,L);  %每层取点
E=exp(r(:,1)*x);    %计算不同h下的场
P=zeros(Lh,1);step=abs(x(end-1));   %dx
for m=1:Lh      %不同h下的P
    P(m)=P(m)+E(m,:)*E(m,:)'*step;
end    
%单独计算第1层的P
x=linspace(0,0.1e-6,L); %每层取点
E=diag(A(:,1))*cos(r(:,2)*x)+diag(B(:,1))*sin(r(:,2)*x);
step=x(2);  %dx
for m=1:Lh     %不同h下的P
    P(m)=P(m)+E(m,:)*E(m,:)'*step;
end    
%单独计算第N+3层的P
x=linspace(0,1e-4,L); %每层取点
E=diag(A(:,end))*exp(-r(:,end)*x);
step=x(2);  %dx
for m=1:Lh      %不同h下的P
    P(m)=P(m)+E(m,:)*E(m,:)'*step;
end
%单独计算第2层的P
for m=1:Lh
    x=linspace(0,b(m,2),L); %每层取点
    E=A(m,2)*cos(r(m,3)*x)+B(m,2)*sin(r(m,3)*x);
    step=x(2);  %dx
    P(m)=P(m)+E*E'*step;    %不同h下的P
end
km=zeros(Lh,1);
h_temp=(-h/N/2)';
for k=4:N+3
    for m=1:Lh
        h_temp(m)=h_temp(m)+h(m)/N;  %每一层的中间线距光栅顶的高度
        d=period*h_temp(m)/h(m);    %计算对应层的w
        Aq=sin(pi*2*d/period)*(3.29^2-3.45^2)/2/pi; %计算Aq
        x=linspace(0,b(m,k-1),L); %每层取点
        if NEFF(m)<n(k) %计算场
            E=A(m,k-1)*cos(r(m,k)*x)+B(m,k-1)*sin(r(m,k)*x);
        else
            E=A(m,k-1)*cosh(r(m,k)*x)+B(m,k-1)*sinh(r(m,k)*x);
        end
        step=x(2);    %dx
        P(m)=P(m)+E*E'*step;    %累加不同h下的P
        km(m)=km(m)+Aq*E*E'*step;   %累加不同h下的km
    end
end
    
km=k0^2/2/beta0*km./P;  %计算耦合系数km
figure;
plot(h,NEFF),xlabel('h/m'),ylabel('n_e_f_f'),%title('n_e_f_f—h曲线');  
figure;
plot(h,km),xlabel('h/m'),ylabel('k_a_b'),%title('k_a_b—h曲线');