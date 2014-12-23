%% clear close
clc, clear all, close all
hold on
axis equal
%% plot rectangle
rectangle('Position',[0.0,0.0,4,1]);

%% draw circle
radius = 0.15;
x0 = 1.5;
y0 = 0.5;
angle = linspace(0,2*pi,50);
x = radius*cos(angle)+x0;
y = radius*sin(angle)+y0;

plot(x,y,'b')




%% plot vector field
% plot upper part
phi1 = pi/2 -acos(1/8);
phi2 = pi/2 -acos(6/8);
phi = linspace(phi1,phi2,10);

x1 = radius*cos(phi)+x0;
y1 = radius*sin(phi)+y0;

%g = @(x,y)  1/(phi2/2-phi1/2)^2*[(x-x0).*(atan2(y-y0,x-x0)-phi1).*(phi2-atan2(y-y0,x-x0)); ...
%                                 (y-y0).*(atan2(y-y0,x-x0)-phi1).*(phi2-atan2(y-y0,x-x0))];
g = @(x,y)   1/(phi2/2-phi1/2)^2*[(x-x0).*( atan((y-y0)./(x-x0))-phi1 ).*( phi2-atan((y-y0)./(x-x0)) ); ...
                                  (y-y0).*( atan((y-y0)./(x-x0))-phi1 ).*( phi2-atan((y-y0)./(x-x0)) )];
gxy = g(x1,y1);

%compute norm of each vector first and last have to be zero and in the middle
%the length should be the radius
gxyuppernorm = sqrt(sum(gxy.^2,1));
quiver(x1,y1,gxy(1,:),gxy(2,:),0,'r')

% plot lower part
phi1 = -pi/2 +acos(1/8);
phi2 = -pi/2 +acos(6/8);
phi = linspace(phi1,phi2,10);

x1 = radius*cos(phi)+x0;
y1 = radius*sin(phi)+y0;

%g = @(x,y)  1/(phi2/2-phi1/2)^2*[(x-x0).*(atan2(y-y0,x-x0)-phi1).*(phi2-atan2(y-y0,x-x0)); ...
%                                 (y-y0).*(atan2(y-y0,x-x0)-phi1).*(phi2-atan2(y-y0,x-x0))];
g = @(x,y)   1/(phi2/2-phi1/2)^2*[(x-x0).*( atan((y-y0)./(x-x0))-phi1 ).*( phi2-atan((y-y0)./(x-x0)) ); ...
                                  (y-y0).*( atan((y-y0)./(x-x0))-phi1 ).*( phi2-atan((y-y0)./(x-x0)) )];
gxy = g(x1,y1);

%compute norm of each vector first and last have to be zero and in the middle
%the length should be the radius
gxylowernorm = sqrt(sum(gxy.^2,1));
quiver(x1,y1,gxy(1,:),gxy(2,:),0,'g')









