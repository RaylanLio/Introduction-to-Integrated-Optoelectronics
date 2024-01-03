clear all clc;

N=20;           %����դ����ΪN��
n=zeros(1,N+4); %n(1)->n0  
neff=linspace(3.291,3.449,1000);
L=length(neff); 
gamma=zeros(L,N+4); %gamma(1)->gamma0  
T=zeros(L,N+3);     %T(1)->T0
A=T;        %ϵ��
B=T;
d=zeros(1,N+2);     %��¼ÿһ��ĺ��
period=250e-9;      %����
wavelength=880e-9;
k0=2*pi/wavelength; 
d(1)=0.1e-6;        %d(1)->d1
NEFF=[];
r=[];
b=[];

n(1)=3.29;      %0��������
n(2)=3.59;      %1��������
n(3)=3.45;      %2��������
n(end)=3.29;    %N+3��������
step_h=(0.1e-6-0.1e-6/1e6)/50;  %ɨ���h�Ĳ���
A_temp=[];
B_temp=[];
h=0.1e-6/1e6:step_h:0.1e-6;     %ɨ���h�ķ�Χ0~0.1e-6
% % disp(['����ʱ��ʼ��']);
for p=1:length(h)
% %     disp(['����',num2str(length(h)-p),'����ɣ�']);
    d(2)=0.15e-6-h(p);    %����2�����
    h0=h(p)/N;  %��դ����֣�ÿ����һ��
    d(3:end)=h0;       
    b=[b;d];        %����ͬh�µ�ÿ���ȴ���b��
    %����ƽ��������
    h_temp=h0*[1:20]-h0/2; %ÿһ����м��߾��դ���ĸ߶�
    x_mid=period*h_temp/h(p);   %�м�
    x2=period-x_mid;
    n(4:N+3)=sqrt((x_mid*3.29^2+x2*3.45^2)/period);%ÿһ���ƽ��������
    err=[];
    for k=1:L
        gamma(k,:)=k0*sqrt(abs(neff(k)^2-n.^2));
        T(k,:)=gamma(k,1:end-1)./gamma(k,2:end);
        A(k,1)=1;         %����A1
        B(k,1)=T(k,1);      %����B1
        for m=2:N+2     %����2~N+2����ϵ��
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
        if neff(k)<n(N+3)   %����N+3����ϵ��
            A(k,N+3)=A(k,N+2)*cos(gamma(k,N+3)*d(N+2))+...
                    B(k,N+2)*sin(gamma(k,N+3)*d(N+2));
        else
            A(k,N+3)=A(k,N+2)*cosh(gamma(k,N+3)*d(N+2))+...
                    B(k,N+2)*sinh(gamma(k,N+3)*d(N+2));
        end
        if neff(k)<n(N+2)   %���㲻ͬneff�µı�������ֵ
        err_temp=tan(gamma(k,N+3)*h0)-(A(k,end-1)+...
            B(k,end-1)*T(k,end))/(A(k,end-1)*T(k,end)-B(k,end-1));
        else
            err_temp=tanh(gamma(k,N+3)*h0)+(A(k,end-1)+...
            B(k,end-1)*T(k,end))/(A(k,end-1)*T(k,end)+B(k,end-1));
        end
        err=[err err_temp];
    end
    err=abs(err);
    flag=find(err==min(err));   %Ѱ�Ҷ�Ӧ�������̽��neff����
    A_temp=[A_temp;A(flag,:)];  %����Ӧ�ĳ�ϵ��������
    B_temp=[B_temp;B(flag,:)];
    r=[r;gamma(flag,:)];        %����Ӧ��gamma������
    NEFF=[NEFF neff(flag)];     
end
figure;
plot(h,NEFF),xlabel('h/m'),ylabel('n_e_f_f')
%title('n_e_f_f��h����');  

A=A_temp;
B=B_temp(:,1:end-1);
save('data.mat','A','B','r','b','h','n','NEFF','N');
clear all,close all,clc;
load data.mat;
period=250e-9;
beta0=2*pi/period;  
k0=2*pi/880e-9;
L=1001;              %ÿһ���ȡ����
Lh=length(h);       %ѡȡ��h�ĵ���
%���������0���P
x=linspace(-1e-4,0,L);  %ÿ��ȡ��
E=exp(r(:,1)*x);    %���㲻ͬh�µĳ�
P=zeros(Lh,1);step=abs(x(end-1));   %dx
for m=1:Lh      %��ͬh�µ�P
    P(m)=P(m)+E(m,:)*E(m,:)'*step;
end    
%���������1���P
x=linspace(0,0.1e-6,L); %ÿ��ȡ��
E=diag(A(:,1))*cos(r(:,2)*x)+diag(B(:,1))*sin(r(:,2)*x);
step=x(2);  %dx
for m=1:Lh     %��ͬh�µ�P
    P(m)=P(m)+E(m,:)*E(m,:)'*step;
end    
%���������N+3���P
x=linspace(0,1e-4,L); %ÿ��ȡ��
E=diag(A(:,end))*exp(-r(:,end)*x);
step=x(2);  %dx
for m=1:Lh      %��ͬh�µ�P
    P(m)=P(m)+E(m,:)*E(m,:)'*step;
end
%���������2���P
for m=1:Lh
    x=linspace(0,b(m,2),L); %ÿ��ȡ��
    E=A(m,2)*cos(r(m,3)*x)+B(m,2)*sin(r(m,3)*x);
    step=x(2);  %dx
    P(m)=P(m)+E*E'*step;    %��ͬh�µ�P
end
km=zeros(Lh,1);
h_temp=(-h/N/2)';
for k=4:N+3
    for m=1:Lh
        h_temp(m)=h_temp(m)+h(m)/N;  %ÿһ����м��߾��դ���ĸ߶�
        d=period*h_temp(m)/h(m);    %�����Ӧ���w
        Aq=sin(pi*2*d/period)*(3.29^2-3.45^2)/2/pi; %����Aq
        x=linspace(0,b(m,k-1),L); %ÿ��ȡ��
        if NEFF(m)<n(k) %���㳡
            E=A(m,k-1)*cos(r(m,k)*x)+B(m,k-1)*sin(r(m,k)*x);
        else
            E=A(m,k-1)*cosh(r(m,k)*x)+B(m,k-1)*sinh(r(m,k)*x);
        end
        step=x(2);    %dx
        P(m)=P(m)+E*E'*step;    %�ۼӲ�ͬh�µ�P
        km(m)=km(m)+Aq*E*E'*step;   %�ۼӲ�ͬh�µ�km
    end
end
    
km=k0^2/2/beta0*km./P;  %�������ϵ��km
figure;
plot(h,NEFF),xlabel('h/m'),ylabel('n_e_f_f'),%title('n_e_f_f��h����');  
figure;
plot(h,km),xlabel('h/m'),ylabel('k_a_b'),%title('k_a_b��h����');