y1=0.9;y2=0.8;y3=1.

n=40000;

figure (1)

figure (1)
for i=1:40
    if i==1
        plot(chain(1:end,i,7),'m-','LineWidth',1); 
    hold on
    elseif i==12
            plot(chain(1:end,i,7),'g-','LineWidth',1); 
    hold on
    elseif i==25
            plot(chain(1:end,i,7),'y-','LineWidth',1); 
    hold on
    elseif i==38
            plot(chain(1:end,i,7),'k-','LineWidth',1); 
    hold on
    else 
        plot(chain(1:end,i,7),'b-','LineWidth',1); 
    end
%     plot(n,y(i),'x','markersize',10,'LineWidth',2);
end
plot([0,n],[1,1],'r--','markersize',15,'LineWidth',2);
% hold on
% plot(n,y2,'xr','markersize',15,'LineWidth',2);
% hold on
% plot(n,y3,'xr','markersize',15,'LineWidth',2);
% hold on
xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
set(gca,'ylim',[0.7 1.2],'ytick',[0.5:0.1:1.2]);

figure(2)
plot(chain(1:end,1,7),'m-','LineWidth',1); 
    hold on
plot(chain(1:end,12,7),'g-','LineWidth',1); 
    hold on
plot(chain(1:end,25,7),'y-','LineWidth',1); 
    hold on
plot(chain(1:end,38,7),'k-','LineWidth',1); 
hold on
plot([0,n],[0.9,0.9],'r--','markersize',15,'LineWidth',2);
plot([0,n],[0.8,0.8],'r--','markersize',15,'LineWidth',2);
legend('\fontsize{15}\bf\theta_1','\fontsize{15}\bf\theta_1_2'...
        ,'\fontsize{15}\bf\theta_2_5','\fontsize{15}\bf\theta_3_8')
xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
set(gca,'ylim',[0.7 1.2],'ytick',[0.5:0.1:1.2])
% legend('\fontsize{15}\bf\theta_1','\fontsize{15}\bf\theta_2','\fontsize{15}\bf\theta_3','\fontsize{15}\bf\theta_4'...
%     ,'\fontsize{15}\bf\theta_5','\fontsize{15}\bf\theta_6','\fontsize{15}\bf\theta_7'...
%     ,'\fontsize{15}\bf\fontsize{15}\bf\theta_8','\fontsize{15}\bf\theta_9','\fontsize{15}\bf\theta_1_0');
% set(legend,'box','off')
Mean=mean(chain(0.7*n:end,:,:));
Std=std(chain(0.7*n:end,:,:));

ParSet = genparset(chain);
Pars = ParSet ( floor ( 0.75 * size(ParSet,1) ) : size(ParSet,1), 1 : 42);
[N1,X1]=density(chain(:,18),[]);
N1=N1./sum(N1); 
figure (4)
plot(X1,N1,'r-','linewidth',1.5)

figure (3)
conver=output.R_stat;
x=conver(2:end,1);
% plot(x,conver(2:end,2),'b-','LineWidth',1); 
% hold on
for i=1:40
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
pos=20000;
plot([pos pos],[1.2 1.65],'r','LineWidth',1.1)
xx=20000;yy=2;
h = text(xx,yy,'1.2');
set(h,'horizontalalignment','center','fontsize',17);
h = plot(xx,yy,'ok');set(h,'markersize',34);

M=Mean(:,:,7)
S=Std(:,:,7)