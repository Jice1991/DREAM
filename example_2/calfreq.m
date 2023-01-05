function F11=calfreq(thE)
global H E L P B
Ns=20; %测试次数
nn=6; %6阶频率
npar=10;
H=0.6;E0=3.3e10;L=10;P=2500;B=0.4;
mu=0.1;
for i=1:npar
   true(:,i)=normrnd(thE(i),mu,Ns,1);
end
A=B*H;
I=B*H*H*H/12;
LL=ones(1,npar);    
LL(:)=L/npar;
% 求解整体的刚度矩阵和质量矩阵；
K=zeros(2*(npar+1));
M=zeros(2*(npar+1));
ModifyM=zeros(2*(npar+1)-2);
ModifyK=zeros(2*(npar+1)-2);
for j=1:Ns
 for R=1:npar
    E=true(j,R)*E0;
    k=Beam1D2Node_Stiffness(E,I,LL(R));
    K=Beam1D2Node_Assembly(K,k,R,R+1);
    m=Beam1D2Node_Mass(P,A,LL(R));
    M=Beam1D2Node_MS_Assembly(M,m,R,R+1);
end
% 采用划行划列法对整体质量矩阵和刚度矩阵进行约束处理；
for i=1:(2*(npar+1)-3)
    for n=1:(2*(npar+1)-3)
        ModifyM(i,n)=M(i+1,n+1);
        ModifyK(i,n)=K(i+1,n+1);
    end
end
ModifyM(2*(npar+1)-2,1:2*(npar+1)-3)=M(2*(npar+1),2:(2*(npar+1)-2));
ModifyM(1:2*(npar+1)-3,2*(npar+1)-2)=M(2:(2*(npar+1)-2),2*(npar+1));
ModifyK(1:2*(npar+1)-3,2*(npar+1)-2)=K(2:(2*(npar+1)-2),2*(npar+1));
ModifyK(2*(npar+1)-2,1:2*(npar+1)-3)=K(2*(npar+1),2:(2*(npar+1)-2));
ModifyM(2*(npar+1)-2,2*(npar+1)-2)=M(2*(npar+1),2*(npar+1));
ModifyK(2*(npar+1)-2,2*(npar+1)-2)=K(2*(npar+1),2*(npar+1));
% 求解简支梁的固有频率和模态振型；
[~,D]=eig(ModifyK,ModifyM);
F=diag(sqrt(D))./(2*pi);
freqtrue(:,j)=F(1:nn);
end
F11=reshape(freqtrue,Ns*nn,1);

