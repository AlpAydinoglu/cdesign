clear all
clc
close all

%data of the experiment (x and \lam values)
load('data.mat')

%parameters
g = 9.81;
mp = 0.1;
mc = 1;
len = 0.5;
alpha = pi/2;
alpha2 = pi/2;
r1 = 0;  
r2 = 0;
d1 = 0.1;
d2 = -0.1;
w1 = d1*sin(alpha2) + len*cos(alpha2) - r1*cos(alpha2);
w2 = -d2*sin(alpha) + len*cos(alpha) - r2*cos(alpha);
Fc = 0.1*eye(2);

%conversion if needed
xhr = xh - [0; pi; 0; 0];
%gain matrices
KK = [3.6931  -46.7030    3.3885   -5.7147];
LL = [-13.9797   13.9767];
j=0;

step = 1;
scale = 1 / ((1/100)*step);

for i = 1:step:length(xh)
subplot(2,2,1)
cla
%SET AXIS
axis([-0.5 0.5 0 0.8])
set(gca,'visible','off')
%RIGHT WALL BACKGROUND
rectangle('Position',[0.1 0.4 0.25 0.2],'FaceColor',[0.5843 0.8157 0.9882],'EdgeColor','none')
hold on
%LEFT WALL BACKGROUND
rectangle('Position',[-0.35 0.4 0.25 0.2],'FaceColor',[0.5843 0.8157 0.9882],'EdgeColor','none')
hold on
%CART
rectangle('Position',[xh(1,i)-(0.1/2) 0 0.1 0.1],'FaceColor',[0 .5 .5])
hold on
%POLE COORDINATES
xpend = xh(1,i)-len*sin(xh(2,i)-pi); %pendulum x coordinate
ypend = len*cos(xh(2,i)-pi); %pendulum y coordinate
%LEFT SPRING
if xpend <= d2
    ne = 4; a = 0.2; ro = 0.05;
    [xs, ys] = spring(-0.35, 0.5, xpend, 0.5,ne,a,ro);
    plot(xs,ys,'LineWidth',2,'Color',[0 0 0]);
    hold on
    %THE WALL THAT TOUCHES THE POLE
    plot([xpend xpend], [0.45 0.55],'LineWidth',3,'Color',[0 0 0])
    hold on
else
    ne = 4; a = 0.2; ro = 0.05;
    [xs, ys] = spring(-0.35, 0.5, d2, 0.5,ne,a,ro);
    plot(xs,ys,'LineWidth',2,'Color',[0 0 0]);
    hold on
    %THE WALL THAT TOUCHES THE POLE
    plot([-0.1 -0.1], [0.45 0.55],'LineWidth',3,'Color',[0 0 0])
    hold on
end
%RIGHT SPRING
if xpend >= d1
    ne = 4; a = 0.2; ro = 0.05;
    [xs, ys] = spring(0.35, 0.5, xpend, 0.5,ne,a,ro);
    plot(xs,ys,'LineWidth',2,'Color',[0 0 0]);
    hold on
    %THE WALL THAT TOUCHES THE POLE
    plot([xpend xpend], [0.45 0.55],'LineWidth',3,'Color',[0 0 0])
    hold on
else
    ne = 4; a = 0.2; ro = 0.05;
    [xs, ys] = spring(0.35, 0.5, d1, 0.5,ne,a,ro);
    plot(xs,ys,'LineWidth',2,'Color',[0 0 0]);
    hold on
    %THE WALL THAT TOUCHES THE POLE
    plot([0.1 0.1], [0.45 0.55],'LineWidth',3,'Color',[0 0 0])
    hold on
end
%WALLS AT THE END OF BACKGROUND
plot([-0.35 -0.35], [0.4 0.6],'LineWidth',5,'Color',[0 0 0])
hold on
plot([0.35 0.35], [0.4 0.6],'LineWidth',5,'Color',[0 0 0])
hold on
%ROD
plot([xh(1,i) xh(1,i)-len*sin(xh(2,i)-pi)], [0.05 len*cos(xh(2,i)-pi)],'LineWidth',2,'Color',[0 0 0])
hold on
%POLE
viscircles([xh(1,i)-len*sin(xh(2,i)-pi) len*cos(xh(2,i)-pi)],0.01,'Color','k','LineWidth',5);
hold on

%calculate \lambda
x1 = xh(1,i); x2 = xh(2,i)-pi; x3 = xh(3,i); x4 = xh(4,i);
p = [x1 - len*sin(x2); len*cos(x2)];
phi1 = d1*sin(alpha) - p(1)*sin(alpha) + p(2)*cos(alpha)-r1*cos(alpha);
phi2 = p(1)*sin(alpha) - d2*sin(alpha) + p(2)*cos(alpha) - r2*cos(alpha);
lambda = pathlcp(Fc,[phi1;phi2]);

u = KK*xhr(:,i)+LL*lambda;
if i > 1
    x1_old = xh(1,i-1); x2_old = xh(2,i-1)-pi;
    p_old = [x1_old - len*sin(x2_old); len*cos(x2_old)];
    phi1_old = d1*sin(alpha) - p_old(1)*sin(alpha) + p_old(2)*cos(alpha)-r1*cos(alpha);
    phi2_old = p_old(1)*sin(alpha) - d2*sin(alpha) + p_old(2)*cos(alpha) - r2*cos(alpha);
    lambda_old = pathlcp(Fc,[phi1_old;phi2_old]);
    u_old = KK*xhr(:,i-1) + LL*lambda_old;
else
    u_old = 0;
end

subplot(2,2,2);
axis([0 6 -20 20])
plot([(j-1)/scale j/scale],[u_old u],'Color','k','LineWidth',2)
xlabel('Time(s)') 
ylabel('u(t)') 
hold on

subplot(2,2,3);
axis([0 6 -0.5 0.5])
plot([(j)/scale (j+1)/scale],[xh(1,i) xh(1,i+1)],'Color','k','LineWidth',2)
xlabel('Time(s)') 
ylabel('x_1(t)') 
hold on

subplot(2,2,4);
axis([0 6 -0.5 0.5])
plot([(j)/scale (j+1)/scale],[xh(2,i)-pi, xh(2,i+1)-pi],'Color','k', 'LineWidth',2)
xlabel('Time(s)') 
ylabel('x_2(t)') 
hold on

j = j + 1;

%GETFRAME
mov(i)=getframe(gcf);
end

mov2 = mov(1:2:end);
%write the video
writerObj = VideoWriter('cartpole_walls', 'MPEG-4');
writerObj.FrameRate = 50;
writerObj.Quality = 100;
open(writerObj);
writeVideo(writerObj, mov2)
close(writerObj);