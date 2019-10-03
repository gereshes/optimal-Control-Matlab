function [statesDot] = dynamicsTemplateForComp(t,states,args,Cont)
%dynamicsTemplateForComp Dynamics for calculating error implementing the
%Clohessy–Wiltshire equations
% 
%Inputs:
%   t - double - current time
%   states - 5 by 1 vector - Current states [x,y,xDot,yDot,m]'
%   args - 3 by 1 - arguments to be passed into the function
%   Cont - 3 by n - Time history of control
%
%Outputs:
%   statesDot - 5 by 1 vector - derivative of state vector
%       
%   Ari Rubinsztejn
%   www.gereshes.com
%   2019.10.02

x=states(1);
y=states(2);
xDot=states(3);
yDot=states(4);
m=states(5);

c1=args(1);
c2=args(2);
n=args(3);

beta=interp1(Cont(3,:),Cont(1,:),t,'spline');
u=interp1(Cont(3,:),Cont(2,:),t,'spline');

xDotDot=(3*n*n*xDot)+(-2*n*yDot)+(u.*sin(beta)./m);
yDotDot=(-2*n*xDot)+(u.*cos(beta)./m);
mDot=-c1*u./c2;
    
statesDot=[xDot;yDot;xDotDot;yDotDot;mDot];
end

