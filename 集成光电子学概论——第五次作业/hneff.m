function [NEFF,km] = hneff(nI,nII,nIII,nIV,tII,tIII,h,Lambda,N,lambda,m)
% Calculate the relationship between h and neff
%   Detailed explanation goes here
%=========每层的实际折射率以及厚度=========
% nI=3.29;tI=inf;
% nII=3.45;tII=0.15;
% nIII=3.59;tIII=0.1;
% nIV=3.29;tIV=inf;
%=========每层的实际折射率=========

%=========光栅参数以及每层的等效折射率=========
% h = 0; %光栅高度
% Lambda = 0.25; %光栅周期
% N = 20; %光栅分层级数
dh = h/N; %光栅层高
n_layer_eff = zeros(N,1);
for i = 1:N
    n_layer_eff(i)=sqrt((i-0.5)/N*nII^2+(1-(i-0.5)/N)*nI^2);
end
n_layer_eff=[nI;n_layer_eff;nII;nIII;nIV];
%=========光栅参数以及每层的等效折射率=========

%=========入射光波长及真空波矢=========
% lambda=0.88; %0.88um 的入射光
k0 = 2*pi/lambda; %真空波矢
%=========入射光波长及真空波矢=========

%=========有效折射率搜寻范围及kappa参数计算=========
neff = [];
dn = 1e-3;
for i=1:N+2
    neff = [neff n_layer_eff(i)+dn:dn:n_layer_eff(i+1)-dn];
end
NUMNEFF = length(neff);
neff = neff';
f = zeros(1,NUMNEFF);
for i=1:NUMNEFF
    i0 = find(n_layer_eff>neff(i),1,"first")-1;
    kappa = k0*sqrt(n_layer_eff.^2-neff(i)^2);
    gamma = [-1i*kappa(1:i0);zeros(N+3-i0,1);-1i*kappa(N+4)];
    kappa = [zeros(i0,1);kappa(i0+1:N+3);0];
    B = zeros(N+3,2);
    B(1,:) = [1 1];
    if i0 < N+1
        if i0 == 1
            B(2,:) = (rotmatrix(kappa(2)*dh)*diag([1,-gamma(1)/kappa(2)])*B(1,:)')';
        else
            B(2,:) = (rotmatrix_h(gamma(2)*dh)*diag([1 -gamma(1)/gamma(2)])*B(1,:)')';
            for j = 2:i0-1
                B(j+1,:) = (rotmatrix_h(gamma(j+1)*dh)*diag([1 gamma(j)/gamma(j+1)])*B(j,:)')';
            end
            B(i0+1,:) = (rotmatrix(kappa(i0+1)*dh)*diag([1 gamma(i0)/kappa(i0+1)])*B(i0,:)')';
        end
        for j=i0+1:N
            B(j+1,:) = (rotmatrix(kappa(j+1)*dh)*diag([1 kappa(j)/kappa(j+1)])*B(j,:)')';
        end
        B(N+2,:) = (rotmatrix(kappa(N+2)*(tII-h))*diag([1 kappa(N+1)/kappa(N+2)])*B(N+1,:)');
        B(N+3,:) = (rotmatrix(kappa(N+3)*tIII)*diag([1 kappa(N+2)/kappa(N+3)])*B(N+2,:)')';
    else
        B(2,:) = (rotmatrix_h(gamma(2)*dh)*diag([1 -gamma(1)/gamma(2)])*B(1,:)')';
        for j = 2:N
            B(j+1,:) = (rotmatrix_h(gamma(j+1)*dh)*diag([1 gamma(j)/gamma(j+1)])*B(j,:))';
        end
        if i0 == N+1
            B(N+2,:) = (rotmatrix(kappa(N+2)*(tII-h))*diag([1 gamma(N+1)/kappa(N+2)])*B(N+1,:)')';
            B(N+3,:) = (rotmatrix(kappa(N+3)*tIII)*diag([1 kappa(N+2)/kappa(N+3)])*B(N+2,:)')';
        else
            B(N+2,:) = (rotmatrix_h(gamma(N+2)*(tII-h))*diag([1 gamma(N+1)/gamma(N+2)])*B(N+1,:)')';
            B(N+3,:) = (rotmatrix(kappa(N+3)*tIII)*diag([1 gamma(N+2)/kappa(N+3)])*B(N+2,:)')';
        end      
    end
    f(i) = gamma(N+4)*B(N+3,1)-kappa(N+3)*B(N+3,2);
    if i > 1
        if f(i)*f(i-1) <= 0
            l1 = -f(i)/(f(i-1)-f(i));
            l0 = f(i-1)/(f(i-1)-f(i));
            NEFF = l1*neff(i-1)+l0*neff(i);

            Nh = 100;
            ddh = dh/Nh;
            xr = linspace(0,dh,Nh+1);
            i0 = find(n_layer_eff(1:N+1)>NEFF,1,"first")-1;
            % i0 = 1;
            kappa = k0*sqrt(n_layer_eff(1:N+1).^2-NEFF^2);
            gamma = [-1i*kappa(1:i0);zeros(N+1-i0,1)];
            kappa = [zeros(i0,1);kappa(i0+1:N+1)];
            
            if i0 == 1
                B(2,:) = (rotmatrix(kappa(2)*dh)*diag([1,-gamma(1)/kappa(2)])*B(1,:)')';
                Ey = B(2,:)*[cos(kappa(2)*xr);sin(kappa(2)*xr)];
                Am = (dh-xr)/h*(nII^2-nI^2).*sinc(pi*m*(dh-xr)/h);
                AEE = Ey*diag(Am)*conj(Ey)'*ddh;
                P = Ey*conj(Ey)'*ddh;
            else
                B(2,:) = (rotmatrix_h(gamma(2)*dh)*diag([1 -gamma(1)/gamma(2)])*B(1,:)')';
                Ey = B(2,:)*[cosh(gamma(2)*xr);sinh(gamma(2)*xr)];
                Am = (dh-xr)/h*(nII^2-nI^2).*sinc(pi*m*(dh-xr)/h);
                AEE = Ey*diag(Am)*conj(Ey)'*ddh;
                P = Ey*conj(Ey)'*ddh;

                for j = 2:i0-1
                    B(j+1,:) = (rotmatrix_h(gamma(j+1)*dh)*diag([1 gamma(j)/gamma(j+1)])*B(j,:)')';
                    Ey = B(j+1,:)*[cosh(gamma(j+1)*xr);sinh(gamma(j+1)*xr)];
                    Am = (j*dh-xr)/h*(nII^2-nI^2).*sinc(pi*m*(j*dh-xr)/h);
                    AEE = AEE+Ey*diag(Am)*conj(Ey)'*ddh;
                    P = P+Ey*conj(Ey)'*ddh;
                end

                if i0 < N+1
                    B(i0+1,:) = (rotmatrix(kappa(i0+1)*dh)*diag([1 gamma(i0)/kappa(i0+1)])*B(i0,:)')';
                    Ey = B(i0+1,:)*[cos(kappa(i0+1)*xr);sin(kappa(i0+1)*xr)];
                    Am = (i0*dh-xr)/h*(nII^2-nI^2).*sinc(pi*m*(i0*dh-xr)/h);
                    AEE = AEE+Ey*diag(Am)*conj(Ey)'*ddh;
                    P = P+Ey*conj(Ey)'*ddh;
                else
                    beta0 = m*pi/Lambda;
                    km = k0^2/(2*beta0*P)*AEE;
                    break;
                end
            end
            for j=i0+1:N
                B(j+1,:) = (rotmatrix(kappa(j+1)*dh)*diag([1 kappa(j)/kappa(j+1)])*B(j,:)')';
                Ey = B(j+1,:)*[cos(kappa(j+1)*xr);sin(kappa(j+1)*xr)];
                Am = (j*dh-xr)/h*(nII^2-nI^2).*sinc(pi*m*(j*dh-xr)/h);
                AEE = AEE+Ey*diag(Am)*conj(Ey)'*ddh;
                P = P+Ey*conj(Ey)'*ddh;
            end
            beta0 = m*pi/Lambda;
            km = k0^2/(2*beta0*P)*AEE;
            break;
        end
    else
        if f(1) == 0
            NEFF = neff(1);
            break;
        end
    end
end
end