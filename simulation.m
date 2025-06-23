clc;clear all;close all;

%设定研究时间范围
T=60;
Ts=0.001;
R=0.45;
Ld=0.00415;
Lq=0.01674;
phir=0.38;
w=1000;

%系统模型
A1=[1-Ts*R/Ld Ts*Lq/Ld*w;-Ts*Ld/Lq*w 1-Ts*R/Lq];
B1=[Ts/Ld 0;0 Ts/Lq];
C1=eye(2);
D1=0;
x1=zeros(2,T+1);
x1(:,1)=[0.8;0.8];

%输入输出维数
p=2;q=2;

%设计测试输入，共161次测试输入
U1Test=[zeros(p*T,1) eye(p*T)];

%得到相应的测试输出
for k=1:p*T+1
    %设置相同初态
    for i=1:T
        u1(:,i)=U1Test(p*i-1:p*i,k);
        x1(:,i+1)=A1*x1(:,i)+B1*u1(:,i);
        Y1Test(q*i-1:q*i,k)=C1*x1(:,i);
    end 
end

    %作为辅助，利用模型计算G1
        %第一对角全为0
        for i=1:T
            G1(q*i-1:q*i,p*i-1:p*i)=zeros(2,2);
        end
        %从第二对角开始计算
        for i=1:T
            for j=1:i-1
                G1(q*i-1:q*i,p*j-1:p*j)=C1*A1^(i-j-1)*B1;
            end
        end
    %利用模型计算L1
        L1=zeros(q*T,2);
        for i=1:T
            L1(q*i-1:q*i,:)=C1*A1^(i-1);
        end
    %利用模型计算Y1Test
    Y1Test_model=G1*U1Test+L1*x1(:,1);

%计算行为相关矩阵
W1=zeros((p+q)*T,p*T);
for i=1:p*T
    W1(:,i)=[U1Test(:,i+1);Y1Test(:,i+1)]-[U1Test(:,1);Y1Test(:,1)];
end
w1off=[U1Test(:,1);Y1Test(:,1)];

%计算可跟踪性集合
    %子空间部分
    WY1=zeros(q*T,p*T);
    for i=1:p*T
        WY1(:,i)=Y1Test(:,i+1)-Y1Test(:,1);
    end
    %偏置部分
    Y1off=Y1Test(:,1);

%设定参考轨迹
sat_limit = 0.6;                % Saturation limit
for i=1:T
    Ydc(q*i-1:q*i,1)=[0.95^(i-1)*sin(pi/6*(i-1));(min(max(sin(pi/8*(i-1)),-sat_limit), sat_limit))'];
    Ydc1(i)=0.95^(i-1)*sin(pi/6*(i-1));
    Ydc2(i)=(min(max(sin(pi/8*(i-1)),-sat_limit), sat_limit))';
end

%判断可跟踪性,track_original为0即可跟踪
g1=pinv(WY1)*(Ydc-Y1off);
track_original=norm(WY1*g1+Y1off-Ydc,2);

%设计互连框架
% K11=[eye(T) zeros(T,T)];
% K12=[zeros(T,T) eye(T)];
% KA1=[zeros(T,T) eye(T)];
% KA2=[zeros(T,T) -eye(T)];
% S1=[zeros(T,T) eye(T)];
% SA=[-eye(T) zeros(T,T)];

K11=[zeros(p*T,p*T) zeros(p*T,p*T)];
K12=[zeros(p*T,p*T) eye(p*T)];
KA1=[eye(p*T) zeros(p*T,p*T)];
KA2=[zeros(p*T,p*T) eye(p*T)];
S1=[eye(p*T) zeros(p*T,p*T)];
SA=[zeros(p*T,p*T) -eye(p*T)];

%设计系统GA
GA=eye(p*T);
WA=orth([eye(p*T);GA]);

%计算相关矩阵
Wsol=null([S1*W1 SA*WA]);
gstar=-pinv([S1*W1 SA*WA])*S1*w1off;
    %判断gstar是否求准
    track_gstar=norm([S1*W1 SA*WA]*gstar+S1*w1off,2);

%计算补偿后的子空间和偏置
WY1comp=[K12*W1 KA2*WA]*Wsol;
Y1offcomp=[K12*W1 KA2*WA]*gstar+K12*w1off;

%判断补偿后的可跟踪性
g1comp=WY1comp\(Ydc-Y1offcomp);
track_compensated=norm(WY1comp*g1comp+Y1offcomp-Ydc,2);


%补偿前系统的迭代学习控制
Ytrack=WY1*g1+Y1off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算数据矩阵K=[k1,k2]
K=Y1Test*pinv([ones(1,p*T+1);U1Test]);
k2=K(:,2:p*T+1);
%预测长度
h=1;

%迭代次数
N=700;
uup=zeros(p*T,N+1);
yup=zeros(p*T,N+1);
%进入迭代
for k=1:N
%求解控制律
    %计算当前迭代输出和误差信息
    for i=1:T
        x1(:,i+1)=A1*x1(:,i)+B1*uup(p*i-1:p*i,k);
        yup(q*i-1:q*i,k)=C1*x1(:,i);
        yup1(i,k)=[1 0]*x1(:,i);
        yup2(i,k)=[0 1]*x1(:,i);
    end
    eup(:,k)=Ydc-yup(:,k);
    %计算Bigk
    for i=1:h
        Bigk(:,:,i)=kron([ones(1,i) zeros(1,h-i)],k2);
    end
    %计算sigmaBigk^TBigk
    sum1=zeros(p*T*h,p*T*h);
    for i=1:h
        sum1=sum1+Bigk(:,:,i)'*Bigk(:,:,i);
    end
    %计算M
    M=sum1+eye(p*T*h);

    %计算sigmaBigk^T
    sum2=zeros(p*T*h,p*T);
    for i=1:h
        sum2=sum2+Bigk(:,:,i)';
    end
    
    %计算N
        N1=sum2;
    %计算输入更新
        DelBigU=M\N1*eup(:,k);
        delu=DelBigU(1:p*T,:);
        uup(:,k+1)=uup(:,k)+delu;
        Enorm(k)=norm(eup(:,k),2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%补偿后系统的迭代学习控制
Ytrackcomp=WY1comp*g1comp+Y1offcomp;
for i=1:T
    Ytrackcomp1(i)=Ytrackcomp(q*i-1);
    Ytrackcomp2(i)=Ytrackcomp(q*i);
end


figure(1)
plot(0:T-1,Ydc1,'-b','LineWidth',4.5);
hold on
plot(0:T-1,yup1(:,N),'-o','LineWidth',4.5);
grid on
plot(0:T-1,Ytrackcomp1,'-*','LineWidth',4.5);
hold on
h=legend('$$y_{1,d}(n)$$','$${y_{1,500}(n)}$$','$${y_{1,500}^{\rm comp}(n)}$$')
set(h,'Interpreter','latex')
hold on
set(gca,'FontName','Times New Roman','FontSize',24,'fontname','Times')
xlabel('Sampling points','Interpreter','latex')
ylabel('Currents of direct axis','Interpreter','latex')

figure(2)
plot(0:T-1,Ydc2,'-b','LineWidth',4.5);
hold on
plot(0:T-1,yup2(:,N),'-o','LineWidth',4.5);
grid on
plot(0:T-1,Ytrackcomp2,'-*','LineWidth',4.5);
hold on
h=legend('$$y_{2,d}(n)$$','$${y_{2,500}(n)}$$','$${y_{2,500}^{\rm comp}(n)}$$')
set(h,'Interpreter','latex')
hold on
set(gca,'FontName','Times New Roman','FontSize',24,'fontname','Times')
xlabel('Sampling points','Interpreter','latex')
ylabel('Currents of quadrature axis','Interpreter','latex')

