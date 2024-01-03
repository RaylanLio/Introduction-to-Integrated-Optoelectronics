clc;
close all;
clear all;

%=========每层的实际折射率以及厚度=========
nI=3.29;tI=inf;
nII=3.45;tII=0.15;
nIII=3.59;tIII=0.1;
nIV=3.29;tIV=inf;
%=========每层的实际折射率=========

%=========光栅参数=========
N1 = 1000;
h = linspace(0,0.1,N1+1); %光栅高度
h = h';
Lambda = 0.25; %光栅周期
N = 20; %光栅分层级数
%=========光栅参数=========

%=========入射光波长=========
lambda=0.88; %0.88um 的入射光
%=========入射光波长=========

NEFF = zeros(N1+1,1);
km = zeros(N1+1,1);
m = 2;
for i=1:N1+1
    [NEFF(i),km(i)]=hneff(nI,nII,nIII,nIV,tII,tIII,h(i),Lambda,N,lambda,m);
    if isnan(km(i))
        km(i) = 0;
    end
end

figure('Name',"h-neff")
plot(h,NEFF);
title("h~n_{eff}");
xlabel("h (\mu m)");
ylabel("n_{eff} (-)");

figure('Name',"h-km")
plot(h,abs(km));
title("h~\kappa_{m}");
xlabel("h (\mu m)");
ylabel("\kappa_{m} (\mu^{-1})");


