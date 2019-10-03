function [cost] = costFunctionU2(XU,shapes)
%costFunctionMf: Cost Function to Minimize the Thrust Squared in Mayer form
%using Simpsons rule
% 
%Assuming thrust is the second entry in the U vector
%
%Inputs:
%   XU - 1 x n vector - vector of the states and actions
%   shapes - 1 x 4 Vector - vector of the X and U shapes
%
%Outputs:
%   cost - double - Cost of the inputs
%
%   Ari Rubinsztejn
%   www.gereshes.com
%   2019.10.02

rx=shapes(1);cx=shapes(2);ru=shapes(3);cu=shapes(4);
x=XU(1:(rx*cx));
u=XU(((rx*cx)+1):end);
X=reshape(x,rx,cx);
u=reshape(u,ru,cu);

T=X(end,:);
hk=diff(T);

xLow=u(2,1:end-1);
xHi =u(2,2:end);


xMid=(xLow+xHi)./2;
cost=(xLow.^2)+(4*(xMid.^2))+(xHi.^2);
cost=cost.*hk/6;
cost=sum(cost);

end

