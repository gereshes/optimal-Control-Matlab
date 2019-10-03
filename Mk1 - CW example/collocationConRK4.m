function [C,Ceq] = collocationConRK4(XU,fun,args,shapes,setConds)
%collocationConRK4 Hermite Simpson Collocation
% 
%Inputs:
%   XU - m by 1 - combined vector of states and controls
%   fun - function handle - function handle of dynamics function
%   args - n by 1 - Vector of arguments for dynamics
%   shapes - 1 by 4 - vector of X and U shapes
%   setConds - n by 2 vector - vector of initial and final conditions
%
%Outputs:
%   C - 6 by 1 vector - nonlinear constraints
%   C - 6 by 1 vector - nonlinear constraints
%       
%
%   Ari Rubinsztejn
%   www.gereshes.com
%   2019.10.02
 
% Unpack the states and controls

rx=shapes(1);cx=shapes(2);ru=shapes(3);cu=shapes(4);
x=XU(1:(rx*cx));
u=XU(((rx*cx)+1):end);
X=reshape(x,rx,cx);
T=X(end,:);

%Set the time bounds 
T(1)=setConds(1,end);
%T(end)=setConds(2,end);

%Calcualte the time step and force even time steps
hk=diff(T);
tHi=T(3:end);
tMid=T(2:end-1);
tLow=T(1:end-2);
tErr=((tMid-tLow)./(tHi-tLow))-.5;
tErr=[tErr,((T(end)-T(end-1))./(T(end)-T(end-2)))-.5];

%set the hard initial and final steps
X(1:5,1)=setConds(1,1:5)';
X(1:4,end)=setConds(2,1:4)';

U=reshape(u,ru,cu);
uLow=U(:,1:end-1);
uHi =U(:,2:end);
uMid=(uLow+uHi)./2;

xLow=X(1:end-1,1:end-1);
xHi =X(1:end-1,2:end);

%RK4 step forward
k1 = feval(fun,xLow,uLow,args);
k2 = feval(fun,xLow + 0.5.*hk.*k1,uMid,args);
k3 = feval(fun,xLow + 0.5.*hk.*k2,uMid,args);
k4 = feval(fun,xLow + hk.*k3,uHi,args);
xNew = xLow + (hk./6).*(k1 + 2*k2 + 2*k3 + k4);


%Packing up the collocation inequality constraints
err= xNew-xHi;
err=[err;tErr];
[n,m] = size(err);
Ceq = (reshape(err,n*m,1));

C = [];


end

