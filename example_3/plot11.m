y1=0.7;y2=0.8;y3=0.9;y4=1.2;

n=15000;
% save('15pars_10cov8damage_2noise','chain','output')
load 15pars_10cov8damage_2noise

figure (1)
npar = 15;
cm = jet(8);
% indx=[1 4 6 10 12 15];
indx=[1 3 6 8 10 11 12 15];
index=1:1:npar;
index(:,indx)=[];
% for i=1:npar
%     set(gca,'colororder',cm)
%     plot(chain(1:end,indx,1),'-','LineWidth',1);  
%     hold on
% end
% hold on
plot(chain(1:end,index,1),'b-','LineWidth',1);
hold on
plot([0,n],[1,1],'r--','markersize',15,'LineWidth',2);
xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
% hold on
% plot(n,y1,'xr','markersize',15,'LineWidth',2);
% hold on
% plot(n,y2,'xr','markersize',15,'LineWidth',2);
% hold on
% plot(n,y3,'xr','markersize',15,'LineWidth',2);
% hold on
% plot(n,y4,'xr','markersize',15,'LineWidth',2);
set(gca,'ylim',[0.5 1.5],'ytick',[0.5:0.2:1.5]);

figure(5)
for i=1:npar
    set(gca,'colororder',cm)
    plot(chain(1:end,indx,1),'-','LineWidth',1);  
    hold on
end
hold on
plot([0,n],[0.9,0.9],'r--','markersize',15,'LineWidth',2);
hold on
plot([0,n],[0.8,0.8],'r--','markersize',15,'LineWidth',2);
hold on
plot([0,n],[0.7,0.7],'r--','markersize',15,'LineWidth',2);
hold on 
plot([0,n],[1.2,1.2],'r--','markersize',15,'LineWidth',2);
xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
set(gca,'ylim',[0.5 1.5],'ytick',[0.5:0.2:1.5]);
% legend('\fontsize{15}\bf\theta_1','\fontsize{15}\bf\theta_2','\fontsize{15}\bf\theta_3','\fontsize{15}\bf\theta_4'...
%     ,'\fontsize{15}\bf\theta_5','\fontsize{15}\bf\theta_6','\fontsize{15}\bf\theta_7'...
%     ,'\fontsize{15}\bf\fontsize{15}\bf\theta_8','\fontsize{15}\bf\theta_9','\fontsize{15}\bf\theta_1_0');
% set(legend,'box','off')
Mean=mean(chain(0.7*n:end,:,:));
Std=std(chain(0.7*n:end,:,:));

ParSet = genparset(chain);
Pars = ParSet ( floor ( 0.75 * size(ParSet,1) ) : size(ParSet,1), 1 : 15);
[N1,X1]=density(chain(:,7),[]);
N1=N1./sum(N1); 
figure (2)
plot(X1,N1,'r-','linewidth',1.5)

figure (3)
conver=output.R_stat;
x=conver(2:end,1);
% plot(x,conver(2:end,2),'b-','LineWidth',1); 
% hold on
for i=1:15
    plot(x,conver(2:end,i+1),'b-','LineWidth',1); 
    hold on
end
% convergence evaluation
plot([0 n*10],[1.2 1.2],'--k','markersize',10,'LineWidth',2);
xlabel('Iteration No.','fontsize',18,'fontweight','bold','fontname','Times');
ylabel('R-statistic','fontsize',18,'fontweight','bold','fontname','Times');
title('Convergence evaluation','fontsize',18,'fontweight','bold','fontname','Times');
set(gca,'fontsize',17);
hold on
pos=100000;
plot([pos pos],[1.2 2],'r','LineWidth',1.1)
xx=100000;yy=2.7;
h = text(xx,yy,'1.2');
set(h,'horizontalalignment','center','fontsize',17);
h = plot(xx,yy,'ok');set(h,'markersize',34);

