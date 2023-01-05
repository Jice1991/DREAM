y1=0.9;y2=0.8;y3=1;

n=40000;
% chain1=chain(:,:,1);chain2=chain(:,:,2);chain3=chain(:,:,3);chain4=chain(:,:,4);
% save ('30pars_1cov8','chain')
% load 30pars_1cov8
load Finaltwospan
% Finaltwospan=chain;
% save('Finaltwospan')
figure (1)
for i=1:30
    if i==2
        plot(chain(1:end,i,4),'m-','LineWidth',1); 
    hold on
    elseif i==6
            plot(chain(1:end,i,4),'g-','LineWidth',1); 
    hold on
    elseif i==13
            plot(chain(1:end,i,4),'y-','LineWidth',1); 
    hold on
    elseif i==19
            plot(chain(1:end,i,4),'k-','LineWidth',1); 
    hold on
    elseif i==24
            plot(chain(1:end,i,4),'c-','LineWidth',1); 
    hold on
    elseif i==28
            plot(chain(1:end,i,4),'r-','LineWidth',1); 
    hold on
    else 
        plot(chain(1:end,i,4),'b-','LineWidth',1); 
    end
%     plot(n,y(i),'x','markersize',10,'LineWidth',2);
end
plot([0,n],[1,1],'r--','markersize',15,'LineWidth',2);
xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
set(gca,'ylim',[0.7 1.2],'ytick',[0.7:0.1:1.2]);

figure(2)
plot(chain(1:end,2,4),'m-','LineWidth',1); 
    hold on
plot(chain(1:end,6,4),'g-','LineWidth',1); 
    hold on
plot(chain(1:end,13,4),'y-','LineWidth',1); 
    hold on
plot(chain(1:end,19,4),'k-','LineWidth',1); 
hold on
plot(chain(1:end,24,4),'c-','LineWidth',1); 
    hold on
plot(chain(1:end,28,4),'b-','LineWidth',1); 
    hold on
plot([0,n],[0.9,0.9],'r--','markersize',15,'LineWidth',2);
plot([0,n],[0.8,0.8],'r--','markersize',15,'LineWidth',2);
set(gca,'ylim',[0.7 1.2],'ytick',[0.7:0.1:1.2]);

xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
% legend('\fontsize{15}\bf\theta_1','\fontsize{15}\bf\theta_2','\fontsize{15}\bf\theta_3','\fontsize{15}\bf\theta_4'...
%     ,'\fontsize{15}\bf\theta_5','\fontsize{15}\bf\theta_6','\fontsize{15}\bf\theta_7'...
%     ,'\fontsize{15}\bf\fontsize{15}\bf\theta_8','\fontsize{15}\bf\theta_9','\fontsize{15}\bf\theta_1_0');
% set(legend,'box','off')
Mean=mean(chain(0.8*n:end,:,:));
Std=std(chain(0.8*n:end,:,:));

ParSet = genparset(chain);
Pars = ParSet ( floor ( 0.7 * size(ParSet,1) ) : size(ParSet,1), 1 : 30);
[N1,X1]=density(chain(:,2),[]);
N1=N1./sum(N1); 
figure (4)
plot(X1,N1,'r-','linewidth',1.5)
basevalue = 0.00;
a=area(X1,N1,basevalue);
newcolors = [0.7 0.7 0.7];colororder(newcolors)
a(1).LineWidth = 1.5;

figure (3)
conver=output.R_stat;
x=conver(2:end,1);
for i=1:30
    plot(x,conver(2:end,i+1),'b-','LineWidth',1); 
    hold on
end
% plot(x,conver(2:end,3),'k-','LineWidth',1); 
% hold on
% plot(x,conver(2:end,4),'g-','LineWidth',1); 
% hold on
% plot(x,conver(2:end,5),'m-','LineWidth',1); 
% hold on
% plot(x,conver(2:end,6),'r-','LineWidth',1); 
% hold on
% plot(x,conver(2:end,7),'b-','LineWidth',1); 
% hold on
% plot(x,conver(2:end,8),'k-','LineWidth',1); 
% hold on
% plot(x,conver(2:end,9),'g-','LineWidth',1); 
% hold on
% plot(x,conver(2:end,10),'m-','LineWidth',1); 
% hold on
% plot(x,conver(2:end,11),'r-','LineWidth',1); 
hold on
% convergence evaluation
plot([0 n*10],[1.2 1.2],'--k','markersize',10,'LineWidth',2);
xlabel('Iteration No.','fontsize',18,'fontweight','bold','fontname','Times');
ylabel('R-statistic','fontsize',18,'fontweight','bold','fontname','Times');
title('Convergence evaluation','fontsize',18,'fontweight','bold','fontname','Times');
set(gca,'fontsize',17);
hold on
pos=40000;
plot([pos pos],[1.2 1.5],'r','LineWidth',1.1)
xx=40000;yy=1.7;
h = text(xx,yy,'1.2');
set(h,'horizontalalignment','center','fontsize',17);
h = plot(xx,yy,'ok');set(h,'markersize',34);

