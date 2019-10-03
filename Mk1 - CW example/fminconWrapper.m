function [xSoln,uSoln,tSoln, Fval, ExitFlag, Output] = fminconWrapper(X,U,fun,args,bds,setConds,linear)
%fminconWrapper Wrapper around fmincon
% 
%Inputs:
%   X - n by m - All n states at every one of m discretization points
%   U - p by m - All p controls at every one of m discretization points
%   fun- function handle - Dynamics of the system
%   args - 1 by c - Vector of c arguments passed into the dynamics
%   bds - n by m Cell array - set of bounds for the states and controls
%   setConds - n by 2 vector - vector of initial and final conditions
%   linear - 1 by 4 Cell array - Set of linear equality and inequality
%       constraints
%
%Outputs:
%   xSoln - n by m - Vector of states 
%   uSoln - p by m - Vector of control 
%   tSoln - 1 by m - Vector of times at those discretization points
%   Fval - double - set of bounds for the states and controls
%   ExitFlag - int - vector of initial and final conditions
%
%   Ari Rubinsztejn
%   www.gereshes.com
%   2019.10.02
 
%Unpack the inputs
[rX,cX]=size(X);
[rU,cU]=size(U);
[rLB,cUB]=size(bds{1});
[rLBC,cUBC]=size(bds{3});
x=reshape(X,rX*cX,1);
u=reshape(U,rU*cU,1);
lbdS=reshape(bds{1},rLB*cUB,1);
ubdS=reshape(bds{2},rLB*cUB,1);
lbdC=reshape(bds{3},rLBC*cUBC,1);
ubdC=reshape(bds{4},rLBC*cUBC,1);
shapes=[rX,cX,rU,cU];
XU=[x;u];
lbd=[lbdS;lbdC];
ubd=[ubdS;ubdC];
A=linear{1};
B=linear{2};
Aeq=linear{3};
Beq=linear{4};


%Choose the cost function and nonlinear Constraints
costFun = @(x)costFunctionMf(x,shapes);
nonlCon = @(x)collocationConRK4(x,fun,args,shapes,setConds);



%Setup the problem object
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','OptimalityTolerance',1e-4,'MaxFunctionEvaluations',1*21000);
Problem.objective = costFun;
Problem.x0 = XU;
Problem.Aineq = A;
Problem.bineq = B;
Problem.Aeq = Aeq;
Problem.beq = Beq;
Problem.lb = lbd;
Problem.ub = ubd;
Problem.nonlcon = nonlCon;
Problem.options = options;
Problem.solver = 'fmincon';

%Plug it into fmincon
[xSoln, Fval, ExitFlag, Output] = fmincon(Problem);


%Unpack the fmincon solution
x=xSoln(1:(rX*cX));
u=xSoln(((rX*cX)+1):end);
X=reshape(x,rX,cX);
tSoln=X(end,:);
xSoln=X(1:end-1,:);
uSoln=reshape(u,rU,cU);
end

