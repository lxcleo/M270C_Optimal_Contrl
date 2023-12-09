max_iteration = 10;
t0 = 0;
tf = 10;
for i = 1:max_iteration
   % 1) start with assumed control u and move forward
   [Tx,X] = ode45(@(t,x) stateEq(t,x,u,Tu), [t0 tf], ...
                  initx, options);
   % 2) Move backward to get the trajectory of costates
   x1 = X(:,1); x2 = X(:,2);
   [Tp,P] = ode45(@(t,p) costateEq(t,p,u,Tu,x1,x2,Tx), ...
                  [tf t0], initp, options);
   p1 = P(:,1);
   % Important: costate is stored in reverse order. The dimension of
   % costates may also different from dimension of states
   % Use interploate to make sure x and p is aligned along the time axis
   p1 = interp1(Tp,p1,Tx);
   % Calculate deltaH with x1(t), x2(t), p1(t), p2(t)
   dH = pH(x1,p1,Tx,u,Tu);
   H_Norm = dH'*dH;
   % Calculate the cost function
   J(i,1) = tf*(((x1')*x1 + (x2')*x2)/length(Tx) + ...
                 0.1*(u*u')/length(Tu));
   % if dH/du < epslon, exit
   if H_Norm < eps
       % Display final cost
       J(i,1)
       break;
   else
       % adjust control for next iteration
       u_old = u;
       u = AdjControl(dH,Tx,u_old,Tu,step);
   end
end
% State equations
function dx = stateEq(t,x,u,Tu)
dx = zeros(2,1);
u = interp1(Tu,u,t); % Interploate the control at time t
dx(1) = -2*(x(1) + 0.25) + (x(2) + 0.5)*exp(25*x(1)/(x(1)+2)) ...
        - (x(1) + 0.25).*u;
dx(2) = 0.5 - x(2) -(x(2) + 0.5)*exp(25*x(1)/(x(1)+2));
end
% Costate equations
function dp = costateEq(t,p,u,Tu,x1,x2,xt)
dp = zeros(2,1);
x1 = interp1(xt,x1,t);   % Interploate the state varialbes
x2 = interp1(xt,x2,t);
u = interp1(Tu,u,t);     % Interploate the control
dp(1) = p(1).*(u + exp((25*x1)/(x1 + 2)).*((25*x1)/(x1 + 2)^2 - ...
        25/(x1 + 2))*(x2 + 1/2) + 2) - ...
        2*x1 - p(2).*exp((25*x1)/(x1 + 2))*((25*x1)/(x1 + 2)^2 - ...
        25/(x1 + 2))*(x2 + 1/2);
dp(2) = p(2).*(exp((25*x1)/(x1 + 2)) + 1) - ...
        p(1).*exp((25*x1)/(x1 + 2)) - 2*x2;
end
% Partial derivative of H with respect to u
function dH = pH(x1,p1,tx,u,Tu)
% interploate the control
u = interp1(Tu,u,tx);
R = 0.1;
dH = 2*R*u - p1.*(x1 + 0.25);
end
% Adjust the control
function u_new = AdjControl(pH,tx,u,tu,step)
% interploate dH/du
pH = interp1(tx,pH,tu);
u_new = u - step*pH;
end