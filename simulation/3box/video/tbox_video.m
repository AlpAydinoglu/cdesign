clear all
clc
close all

%data of the experiment (x and \lam values)
load('data.mat')

%parameters
mc = 1;
m1 = 1;
m2 = 1;
m3 = 1;
g = 9.81;
mp = 1.5;
len = 0.5; %0.5 to 0.2
Fc = 0.01*eye(2);

%conversion if needed
xhr = xh - [0; 0; 0; pi; 0; 0; 0; 0];
%gain matrices
KK = [5.2182    0.0000    0.6848 -353.2109    8.9340    0.0000   -3.1900  -55.1658; 6.8662   -0.0000   -6.5382  -51.2013    2.3006   -0.0000   -2.7591  -10.0627];
LL = [-4.0164   -0.2520; -0.2520    6.4276];
j=0;

d = 1;
T = length(xh);
mov = struct('cdata', cell(1,T+1), 'colormap', cell(1,T+1));
xh_h = xh;
%xh = xh_h(1:8,1:10:end);
for i = 1:length(xh)
%for i = 1:5000
    if mod(i,10) == 0
        %i
    end
%f = figure;
%subplot(2,2,1)
cla
%SET AXIS
axis([-3 3 0 0.6])
set(gca,'visible','off')
set(gcf, 'Position',  [100, 100, 1000, 500])
%CART 1
wc1 = 1; %cart 1 width changd from 0.2 to 1
lc1 = 0.05; %cart 1 length
rectangle('Position',[xh(1,i)-(wc1/2)-d 0 wc1 lc1],'FaceColor',[0 .5 .5])
hold on
%POLE COORDINATES
xpend = xh(1,i)-len*sin(xh(4,i)-pi)-d; %pendulum x coordinate
ypend = len*cos(xh(4,i)-pi); %pendulum y coordinate
%ROD
plot([xh(1,i)-d xh(1,i)-d-len*sin(xh(4,i)-pi)], [lc1/2 len*cos(xh(4,i)-pi)],'LineWidth',2,'Color',[0 0 0])
hold on
%POLE
%viscircles([xh(1,i)-d-len*sin(xh(4,i)-pi) len*cos(xh(4,i)-pi)],0.05,'Color','k','LineWidth',3);
plot( xh(1,i)-d-len*sin(xh(4,i)-pi), len*cos(xh(4,i)-pi), '.k', 'MarkerSize',50)
hold on
%CART 2
wc2 = 0.7;
lc2 = 0.05;
rectangle('Position',[xh(2,i)-(wc2/2) 0 wc2 lc2],'FaceColor',[0 .5 .5])
hold on
%CART 3
wc3 = 1;
lc3 = 0.05;
rectangle('Position',[xh(3,i)+d-(wc3/2) 0 wc3 lc3],'FaceColor',[0 .5 .5])
hold on
%RESPECTIVE POINTS
%plot(-d,0,'-x','Color','r','LineWidth',2)
%rectangle('Position',[-d-0.05 0 0.1 0.01],'FaceColor', 'k')
%patch([-d-0.05 -d-0.05 -d+0.05 -d+0.05], [0 0.01 0.01 0], 'red', 'FaceAlpha', 0.2)
plot([-d -d],[0, 0.013], 'LineWidth', 2, 'Color', 'k')
text(-d-0.2,-0.02,'x_1= 0')
hold on
%plot(0,0,'-x','Color','r','LineWidth',2)
%patch([0-0.05 0-0.05 0+0.05 0+0.05], [0 0.01 0.01 0], 'red', 'FaceAlpha', 0.2)
plot([0 0],[0, 0.013], 'LineWidth', 2, 'Color', 'k')
text(-0.2,-0.02,'x_2= 0')
hold on
%plot(d,0,'-x','Color','r','LineWidth',2)
%patch([d-0.05 d-0.05 d+0.05 d+0.05], [0 0.01 0.01 0], 'red', 'FaceAlpha', 0.2)
plot([d d],[0, 0.013], 'LineWidth', 2, 'Color', 'k')
text(d-0.2,-0.02,'x_3= 0')
hold on
%PLOT CONNECTIONS
%FOR BOX1
% plot( [xh(1,i)-d xh(1,i)-d/2], [lc1/2 lc1/2],'Color','r', 'LineWidth', 3);
% hold on
% plot( [xh(1,i)-d/2 xh(1,i)-d/2], [0.25*lc1 0.75*lc1],'Color','r', 'LineWidth', 3);

%FOR BOX2
% plot( [xh(2,i) xh(2,i)+d/2], [0.05 0.05],'Color','k', 'LineWidth', 2);
% hold on
% plot( [xh(2,i) xh(2,i)-d/2], [0.05 0.05],'Color','k', 'LineWidth', 2);
% hold on

%FOR BOX3
% plot( [xh(3,i)+d xh(3,i)+d-d/2], [lc3/2 lc3/2],'Color','r', 'LineWidth', 3);
% hold on
% plot( [xh(3,i)+d-d/2 xh(3,i)+d-d/2], [0.25*lc3 0.75*lc3],'Color','r', 'LineWidth', 3);
% hold on

ne = 7;  ro = 0.01;

%SPRINGS FOR BOX1
a = d/2;
dist1 = xh(2,i) - xh(1,i); %distance between x2 and x1
if dist1 <= 0
    [xs, ys] = spring(xh(1,i)-d/2, lc2/2, xh(2,i)-d/10, lc2/2,ne,a,ro);
    plot(xs,ys,'LineWidth',2,'Color',[0 0 0]);
    hold on
%     plot( [xh(2,i)+d/2 xh(2,i)+d/2], [0.25*lc2 0.75*lc2],'Color','k', 'LineWidth', 3);
%     hold on
else
    [xs, ys] = spring(xh(2,i)-d/2, lc2/2, xh(2,i)-d/10, lc2/2,ne,a,ro);
    plot(xs,ys,'LineWidth',2,'Color',[0 0 0]);
    hold on
%     plot( [xh(2,i)+d/2 xh(2,i)+d/2], [0.25*lc2 0.75*lc2],'Color','k', 'LineWidth', 3);
%     hold on
end

%SPRINGS FOR BOX2
a = d/2;
dist2 = xh(3,i) - xh(2,i); %distance between x3 and x2
if dist2 <= 0
    [xs, ys] = spring(xh(2,i)+d/10, lc2/2, xh(3,i)+d/2, lc2/2,ne,a,ro);
    plot(xs,ys,'LineWidth',2,'Color',[0 0 0]);
    hold on
%     plot( [xh(2,i)+d/2 xh(2,i)+d/2], [0.25*lc2 0.75*lc2],'Color','k', 'LineWidth', 3);
%     hold on
else
    [xs, ys] = spring(xh(2,i)+d/10, lc2/2, xh(2,i)+d/2, lc2/2,ne,a,ro);
    plot(xs,ys,'LineWidth',2,'Color',[0 0 0]);
    hold on
%     plot( [xh(2,i)+d/2 xh(2,i)+d/2], [0.25*lc2 0.75*lc2],'Color','k', 'LineWidth', 3);
%     hold on
end

%calculate \lambda
q1 = xh(1,i); q3 = xh(2,i); q4 = xh(3,i);
phi1 = q3 - q1;
phi2 = q4 - q3;
lambda = pathlcp(Fc,[phi1;phi2]);


u = KK*xhr(:,i)+LL*lambda;
if i > 1
    q1 = xh(1,i-1); q3 = xh(2,i-1); q4 = xh(3,i-1);
    phi1 = q3 - q1;
    phi2 = q4 - q3;
    lambda_old = pathlcp(Fc,[phi1;phi2]);
    u_old = KK*xhr(:,i-1) + LL*lambda_old;
else
    u_old = [0;0];
end
    
% subplot(2,2,3);
% axis([0 50 -60 60])
% plot([(j-1)/100 j/100],[u_old(1) u(1)],'Color','k')
% xlabel('Time(s)') 
% ylabel('u_1(t)') 
% hold on

% subplot(2,2,2);
% %axis([0 20 -20 20])
% plot([(j-1)/100 j/100],[u_old(2) u(2)],'Color','k')
% xlabel('Time(s)') 
% ylabel('u_2(t)-zm') 
% hold on

% subplot(2,2,4);
% axis([0 50 -60 60])
% plot([(j-1)/100 j/100],[u_old(2) u(2)],'Color','k')
% xlabel('Time(s)') 
% ylabel('u_2(t)') 
% hold on

% subplot(2,2,4);
% %axis([0 6 -0.5 0.5])
% plot(j/100,xh(2,i)-pi,'.','Color','k')
% xlabel('Time(s)') 
% ylabel('x_2(t)') 
% hold on

j = j + 1;

%GETFRAME
mov(i)=getframe(gcf);
end

%write the video
mov2 = mov(1:2:2500);
writerObj = VideoWriter('3box_example_rt', 'MPEG-4');
writerObj.FrameRate = 50;
writerObj.Quality = 100;
open(writerObj);
writeVideo(writerObj, mov2)
close(writerObj);