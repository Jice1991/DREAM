
%%
%Two-span continuous beam
clear;clc;
global H E L P B
Ns=20; % measurement No.
nn=8; % frequencies No.
mm=8;% mode shape No.
npar=30;
H=0.25;E0=2.8e7;L=3.8;P=2300;B=0.15; % structural property
%exact theta 
thetaE=ones(1,npar);
% mu=0.05; % adjust mu to vary measurement error
for i=1:npar
thetaE(i)=1+0.*i;
end
% thetaE=[0.7 1 1 1 0.7 1 1 1 1 1 0.8 1 0.8 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
thetaE(2)=0.8;thetaE(6)=0.9;thetaE(13)=0.8;
thetaE(19)=0.8;thetaE(24)=0.9;thetaE(28)=0.9;
% for i=1:npar
%    true(:,i)=normrnd(thetaE(i),mu,Ns,1);
% end
true=thetaE;
A=B*H;
I=B*H*H*H/12;
LL=ones(1,npar);    
LL(:)=L/npar;
% 求解整体的刚度矩阵和质量矩阵；
K=zeros(2*(npar+1));
M=zeros(2*(npar+1));
% ModifyM=zeros(2*(npar+1)-2);
% ModifyK=zeros(2*(npar+1)-2);
for j=1:1 % only measure once
 for R=1:npar
    E=true(j,R)*E0;
    k=Beam1D2Node_Stiffness(E,I,LL(R));
    K=Beam1D2Node_Assembly(K,k,R,R+1);
    m=Beam1D2Node_Mass(P,A,LL(R));
    M=Beam1D2Node_MS_Assembly(M,m,R,R+1);
end
% 采用划行划列法对整体质量矩阵和刚度矩阵进行约束处理；
% for i=1:(2*(npar+1)-3)
%     for n=1:(2*(npar+1)-3)
%         ModifyM(i,n)=M(i+1,n+1);
%         ModifyK(i,n)=K(i+1,n+1);
%     end
% end
% ModifyM(2*(npar+1)-2,1:2*(npar+1)-3)=M(2*(npar+1),2:(2*(npar+1)-2));
% ModifyM(1:2*(npar+1)-3,2*(npar+1)-2)=M(2:(2*(npar+1)-2),2*(npar+1));
% ModifyK(1:2*(npar+1)-3,2*(npar+1)-2)=K(2:(2*(npar+1)-2),2*(npar+1));
% ModifyK(2*(npar+1)-2,1:2*(npar+1)-3)=K(2*(npar+1),2:(2*(npar+1)-2));
% ModifyM(2*(npar+1)-2,2*(npar+1)-2)=M(2*(npar+1),2*(npar+1));
% ModifyK(2*(npar+1)-2,2*(npar+1)-2)=K(2*(npar+1),2*(npar+1));
ModifyM=M;
ModifyK=K;
ModifyM([1,npar+1,2*npar+1],:)=[];
ModifyM(:,[1,npar+1,2*npar+1])=[];
ModifyK([1,npar+1,2*npar+1],:)=[];
ModifyK(:,[1,npar+1,2*npar+1])=[];

% find eigenvector and eigenvalue
[V,D]=eig(ModifyK,ModifyM);
F=diag(sqrt(D))./(2*pi);
% unit length normalization
ModifyV=[V(:,1)./V(1,1) V(:,2)./V(1,2) V(:,3)./V(1,3) V(:,4)./V(1,4) V(:,5)./V(1,5) V(:,6)./V(1,6)...
    V(:,7)./V(1,7) V(:,8)./V(1,8) V(:,9)./V(1,9) V(:,10)./V(1,10)];
%remove rotational DOF
% ModifyV1=ModifyV(2:2:(size(ModifyV,1)-1),:);
No = size(ModifyV,1);
ModifyV1=ModifyV([2:2:(No-1)/2,((No+1)/2+1):2:(No-1)],:);
freqtrue(:,j)=F(1:nn);
modeltrue1(:,:,j)=ModifyV1(:,1:mm);
end
freqmu=freqtrue;modelmu=modeltrue1;
%% add noise
for i=1:Ns
    freqtrue(:,i)=freqmu.*(1+unifrnd(-0.013,0.013,8,1)); %  frequencies-- 1%cov
    modeltrue1(:,:,i)=modelmu.*(1+unifrnd(-0.013,0.013,28,8));%  mode shapes
end

%%
s1=var(freqtrue',1);
modeltrue=remodel(modeltrue1);% normalization
for j=1:mm
    for i=1:Ns
        aa(:,i)=modeltrue(:,j,i);
    end
    mu=mean(aa')';
    mu=repmat(mu,1,Ns);
    kk=aa-mu;
    for i=1:Ns
        kkk(i)=(norm(kk(:,i)))^2;
    end
    s2(j)=sum(kkk)/Ns;
end
% save ('freqtrue_2span','freqtrue')
% save ('modeltrue_2span','modeltrue')
% save ('s1_2span','s1')
% save ('s2_2span','s2')

freqmu=mean(freqtrue,2);
%画图
for jj=1:nn
    subplot(4,3,jj)
    plot(1:Ns,freqtrue(jj,:)','-bo',[0,Ns],[freqmu(jj) freqmu(jj)],'r--','linewidth',2)
end
% legend('测试值','平均值')
