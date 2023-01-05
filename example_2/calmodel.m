function F11=calmodel(thE)
global H E L P B
Ns=20; %���Դ���
mm=3;%3������
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
% �������ĸնȾ������������
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
% ���û��л��з���������������͸նȾ������Լ������
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
% ����֧���Ĺ���Ƶ�ʺ�ģ̬���ͣ�
[V,D]=eig(ModifyK,ModifyM);
% �������������д���
ModifyV=[V(:,1)./V(1,1) V(:,2)./V(1,2) V(:,3)./V(1,3) V(:,4)./V(1,4)];
%ż���д����Ŷ�
ModifyV1=ModifyV(2:2:(size(ModifyV,1)-1),:);
modeltrue1(:,:,j)=ModifyV1(:,1:mm);
end
modeltrue=remodel(modeltrue1);
F11=reshape(reshape(modeltrue,(npar-1)*mm,Ns),(npar-1)*mm*Ns,1);