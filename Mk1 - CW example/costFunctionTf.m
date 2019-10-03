function [cost] = costFunctiontTf(XU,shapes)
%Cost Function to Minimize the Final Time in Meyer form using Simpsons rule
%
%Assuming time is the last entry in the X vector
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
X=reshape(x,rx,cx);



cost=-X(end);
end

