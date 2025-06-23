clc;clear all;close all;

%simulation parameters--unavailable in simulations
T=60;
Ts=0.001;
R=0.45;
Ld=0.00415;
Lq=0.01674;
phir=0.38;
w=1000;

%state-space model
A1=[1-Ts*R/Ld Ts*Lq/Ld*w;-Ts*Ld/Lq*w 1-Ts*R/Lq];
B1=[Ts/Ld 0;0 Ts/Lq];
C1=eye(2);
D1=0;
x1=zeros(2,T+1);
x1(:,1)=[0.8;0.8];

%input-output dimension
p=2;q=2;

%test inputs for 161 times
U1Test=[zeros(p*T,1) eye(p*T)];

%test outputs
for k=1:p*T+1
    for i=1:T
        u1(:,i)=U1Test(p*i-1:p*i,k);
        x1(:,i+1)=A1*x1(:,i)+B1*u1(:,i);
        Y1Test(q*i-1:q*i,k)=C1*x1(:,i);
    end 
end

%     %the matrix G1
%         for i=1:T
%             G1(q*i-1:q*i,p*i-1:p*i)=zeros(2,2);
%         end
%         for i=1:T
%             for j=1:i-1
%                 G1(q*i-1:q*i,p*j-1:p*j)=C1*A1^(i-j-1)*B1;
%             end
%         end
%     %The matrix L1
%         L1=zeros(q*T,2);
%         for i=1:T
%             L1(q*i-1:q*i,:)=C1*A1^(i-1);
%         end


%admissible bahavior-subspace and offset component
W1=zeros((p+q)*T,p*T);
for i=1:p*T
    W1(:,i)=[U1Test(:,i+1);Y1Test(:,i+1)]-[U1Test(:,1);Y1Test(:,1)];
end
w1off=[U1Test(:,1);Y1Test(:,1)];

%original trackability set
    %subspace
    WY1=zeros(q*T,p*T);
    for i=1:p*T
        WY1(:,i)=Y1Test(:,i+1)-Y1Test(:,1);
    end
    %offset
    Y1off=Y1Test(:,1);

%desired reference
sat_limit = 0.6;     % Saturation limit
for i=1:T
    Ydc(q*i-1:q*i,1)=[0.95^(i-1)*sin(pi/6*(i-1));(min(max(sin(pi/8*(i-1)),-sat_limit), sat_limit))'];
    Ydc1(i)=0.95^(i-1)*sin(pi/6*(i-1));
    Ydc2(i)=(min(max(sin(pi/8*(i-1)),-sat_limit), sat_limit))';
end

%verification of trackability, Ydc is trackable when track_original=0
g1=pinv(WY1)*(Ydc-Y1off);
track_original=norm(WY1*g1+Y1off-Ydc,2);

%interconnection-related parameters
K11=[zeros(p*T,p*T) zeros(p*T,p*T)];
K12=[zeros(p*T,p*T) eye(p*T)];
KA1=[eye(p*T) zeros(p*T,p*T)];
KA2=[zeros(p*T,p*T) eye(p*T)];
S1=[eye(p*T) zeros(p*T,p*T)];
SA=[zeros(p*T,p*T) -eye(p*T)];

%auxiliary system SigmaA
GA=eye(p*T);
WA=orth([eye(p*T);GA]);

%relevant matrices
Wsol=null([S1*W1 SA*WA]);
gstar=-pinv([S1*W1 SA*WA])*S1*w1off;
%     track_gstar=norm([S1*W1 SA*WA]*gstar+S1*w1off,2);

%compensated trackability set-subspace and offset
WY1comp=[K12*W1 KA2*WA]*Wsol;
Y1offcomp=[K12*W1 KA2*WA]*gstar+K12*w1off;

%verification of trackability-after compensation, Y_d^c is trackable when
%track_compensated=0
g1comp=WY1comp\(Ydc-Y1offcomp);
track_compensated=norm(WY1comp*g1comp+Y1offcomp-Ydc,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data-based ILC for original system
K=Y1Test*pinv([ones(1,p*T+1);U1Test]);
k2=K(:,2:p*T+1);
%预测长度
h=1;

%iteration number
N=700;
uup=zeros(p*T,N+1);
yup=zeros(p*T,N+1);
for k=1:N
%solution of control law
    for i=1:T
        x1(:,i+1)=A1*x1(:,i)+B1*uup(p*i-1:p*i,k);
        yup(q*i-1:q*i,k)=C1*x1(:,i);
        yup1(i,k)=[1 0]*x1(:,i);
        yup2(i,k)=[0 1]*x1(:,i);
    end
    eup(:,k)=Ydc-yup(:,k);
    for i=1:h
        Bigk(:,:,i)=kron([ones(1,i) zeros(1,h-i)],k2);
    end
    sum1=zeros(p*T*h,p*T*h);
    for i=1:h
        sum1=sum1+Bigk(:,:,i)'*Bigk(:,:,i);
    end
    M=sum1+eye(p*T*h);

    sum2=zeros(p*T*h,p*T);
    for i=1:h
        sum2=sum2+Bigk(:,:,i)';
    end

        N1=sum2;
        DelBigU=M\N1*eup(:,k);
        delu=DelBigU(1:p*T,:);
        uup(:,k+1)=uup(:,k)+delu;
        Enorm(k)=norm(eup(:,k),2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tracking performance after compensation
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

