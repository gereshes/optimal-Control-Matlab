function [statesDot] = dynamicsTemplate(states,control,args)
%dynamicsTemplate Dynamics for collocation implementing the
%Clohessy–Wiltshire equations
% 
%Inputs:
%   states - 5 by m vector - Current states [x,y,xDot,yDot,m]'
%   control - 2 by m - Time history of control
%   args - 3 by 1 - arguments to be passed into the function
%
%Outputs:
%   statesDot - 5 by m vector - derivative of state vector
%       
%   Ari Rubinsztejn
%   www.gereshes.com
%   2019.10.02

x=states(1,:);
y=states(2,:);
xDot=states(3,:);
yDot=states(4,:);
m=states(5,:);

c1=args(1);
c2=args(2);
n=args(3);

beta=control(1,:);
u=control(2,:);

xDotDot=(3*n*n*xDot)+(-2*n*yDot)+(u.*sin(beta)./m);
yDotDot=(-2*n*xDot)+(u.*cos(beta)./m);
mDot=-c1*u./c2;
    
statesDot=[xDot;yDot;xDotDot;yDotDot;mDot];
end

