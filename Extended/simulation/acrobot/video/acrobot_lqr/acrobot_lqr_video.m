clear all
clc
close all

%data of the experiment (x and \lam values)
load('data.mat')

%parameters
len1 = 0.5; %length of the pole 1
len2 = 1; %length of the pole 2
angle_cons = 0.2; %angle in radian
d1 = angle_cons;
d2 = -angle_cons;
xh = x_save;

KK = [1476.33342658263,851.675691160606,548.808564426145,334.434444770683];
LL = [0 0];

counter = 0;
j=0;

for i = 41:20:length(xh)
counter = counter+1;
subplot(2,2,1);
cla
%SET AXIS
axis([-2 2 0 1.6])
set(gca,'visible','off')

p1x = len1*sin(xh(1,i)); %pendulum 1 x coordinate
p1y = len1*cos(xh(1,i)); %pendulum y coordinate

p2x = p1x - len2*sin(xh(2,i));
p2y = p1y + len2*cos(xh(2,i));

%ROD1
plot([0 p1x], [0 p1y],'LineWidth',1,'Color',[0 0 1])
hold on
%ROD2
plot([p1x p2x], [p1y p2y],'LineWidth',1,'Color',[0 0 1])
hold on
%LIMITS
% plot([0 len1*sin(d1)], [0 len1*cos(d1)],'LineWidth',1,'Color',[1 0 0])
% plot([0 len1*sin(d2)], [0 len1*cos(d2)],'LineWidth',1,'Color',[1 0 0])
patch_x1 = [0 len1*sin(d1) 2];
patch_y1 = [0 len1*cos(d1) 0];
patch(patch_x1,patch_y1,'black','FaceAlpha',0.1,'EdgeColor','none')
hold on
patch_x2 = [0 len1*sin(d2) -2];
patch_y2 = [0 len1*cos(d2) 0];
patch(patch_x2,patch_y2,'black','FaceAlpha',0.1,'EdgeColor','none')
%POLES
viscircles([p1x p1y],0.01,'Color','r','LineWidth',5);
hold on
%viscircles([p2x p2y],0.02,'Color','k','LineWidth',5);
%hold on
viscircles([0 0],0.04,'Color','k','LineWidth',5);
hold on

%Calculate lambda
x1 = xh(1,i); x2 = xh(2,i); x3 = xh(3,i); xh4 = xh(4,i);
phi1 = d1 - x1;
phi2 = x1 - d2;
Fc = eye(2);
lambda = pathlcp(Fc,[phi1;phi2]);

u = KK*xh(:,i)+LL*lambda;
if i > 1
    phi1 = d1 - xh(1,i-20);
    phi2 = xh(1,i-20) - d2;
    lambda_old = pathlcp(Fc,[phi1;phi2]);
    u_old = KK*xh(:,i-20) + LL*lambda_old;
else
    u_old = 0;
end

subplot(2,2,2);
axis([0 2.2 -30 2])
plot([(j-1)/50 j/50],[u_old u],'Color','k', 'LineWidth', 2)
xlabel('Time(s)') 
ylabel('u(t)') 
hold on

subplot(2,2,3);
axis([0 2.2 -2 2])
plot(j/50,xh(1,i),'.','Color','k')
xlabel('Time(s)') 
ylabel('x_1(t)') 
hold on

subplot(2,2,4);
axis([0 2.2 -2 2])
plot(j/50,xh(2,i),'.','Color','k')
xlabel('Time(s)') 
ylabel('x_2(t)') 
hold on

j = j+1;
%TEXTS
%txt = sprintf('l_1 = %f', lamh(1,i));
%text(0.15,0.6,txt);
%hold on
%txt2 = sprintf('l_2 = %f', lamh(2,i));
%text(-0.33,0.6,txt2);
%hold on
%KK = [3.6931  -46.7030    3.3885   -5.7147];
%LL = [-13.9797   13.9767];
%a1 = KK*xhr(:,i);
%a2 = LL*lamh(:,i);
%txt3 = sprintf('u_x + u_l = (%f) + (%f) = (%f)', a1,a2,a1+a2);
%text(-0.2,0.9,txt3);
%hold on
%txt4 = sprintf('theta = %f', xh(2,i));
%text(-0.2,0.8,txt4);
%hold on
%txt5 = sprintf('thetadot = %f', xh(4,i));
%text(-0.2,0.75,txt5);

%GETFRAME
mov(counter)=getframe(gcf);
end

%write the video
writerObj = VideoWriter('acrobot_lqr_fail', 'MPEG-4');
writerObj.FrameRate = 50;
writerObj.Quality = 100;
open(writerObj);
writeVideo(writerObj, mov)
close(writerObj);