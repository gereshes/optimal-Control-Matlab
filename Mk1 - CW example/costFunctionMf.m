function [cost] = costFunctionMf(XU,shapes)
%costFunctionMf: Cost Function to Minimize the Final Mass in Mayer Form using Simpsons rule
%
%Assuming mass is the second to last entry in the X vector
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


cost=X(end-1,end);

end

