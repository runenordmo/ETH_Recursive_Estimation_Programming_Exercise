function [posEst,linVelEst,oriEst,driftEst,...
          posVar,linVelVar,oriVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,driftEst,...
%    posVar,linVelVar,oriVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConstants.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2018
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch
%

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % initial state mean
    posEst = [0 0];
    linVelEst = [0 0];
    oriEst = 0;
    driftEst = 0;
    
    % initial state variance
    posVar = [(1/4)*estConst.StartRadiusBound^2 (1/4)*estConst.StartRadiusBound^2];
    linVelVar = [0 0];
    oriVar = (1/3)*estConst.RotationStartBound.^2;
    driftVar = 0;
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar,linVelVar,oriVar,driftVar]);
    % estimator state
    estState.xm = [posEst'; linVelEst'; oriEst; driftEst];
    % time of last update
    estState.tm = tm;
    return;
end

%% Estimator iteration.

% get time since last estimator update
dt = tm - estState.tm;
estState.tm = tm; % update measurement update time

% prior update
% solve differential equation for mean x_prior[k]
xpAndPp0 = [estState.xm; reshape(estState.Pm,[36 1])];
[t,xpAndPp] = ode45(@(t,xpAndPp) updateDiffEq(t,xpAndPp,actuate,estConst),[tm-dt tm], [xpAndPp0]);
xpAndPpEnd = xpAndPp(end,:)';

xp = xpAndPpEnd(1:6);
Pp = reshape(xpAndPpEnd(7:42),[6 6]);

% measurement update

% Set up convenient variables
xa = estConst.pos_radioA(1);
ya = estConst.pos_radioA(2);
xb = estConst.pos_radioB(1);
yb = estConst.pos_radioB(2);
xc = estConst.pos_radioC(1);
yc = estConst.pos_radioC(2);

xDistA = xp(1) - xa;
yDistA = xp(2) - ya;
xDistB = xp(1) - xb;
yDistB = xp(2) - yb;
xDistC = xp(1) - xc;
yDistC = xp(2) - yc;
distA = norm([xp(1);xp(2)]-[xa;ya],2);
distB = norm([xp(1);xp(2)]-[xb;yb],2);
distC = norm([xp(1);xp(2)]-[xc;yc],2);

% Set equal to prior update initially
xm = xp; Pm = Pp;

% Extract measurement
zk = sense';

% Perform measurement update according to availability of measurement 3 
measurement3Available = isfinite(zk(3));
if(measurement3Available)
    Hk = zeros(5,6);
    Hk = [ ...
        xDistA/distA yDistA/distA 0 0 0 0;
        xDistB/distB yDistB/distB 0 0 0 0;
        xDistC/distC yDistC/distC 0 0 0 0;
                   0            0 0 0 1 1;
                   0            0 0 0 1 0;
    ];
    Mk = eye(5);
    hk_xp = [ ...
        distA;
        distB;
        distC;
        xp(5) + xp(6);
        xp(5);
    ];
    R = diag([estConst.DistNoiseA,estConst.DistNoiseB, ...
        estConst.DistNoiseC, estConst.GyroNoise, estConst.CompassNoise]);
else
    zk = [zk(1:2);zk(4:5)];
    Hk = zeros(4,6);
    Hk = [ ...
        xDistA/distA yDistA/distA 0 0 0 0;
        xDistB/distB yDistB/distB 0 0 0 0;
                   0            0 0 0 1 1;
                   0            0 0 0 1 0;
    ];
    Mk = eye(4);
    hk_xp = [ ...
        distA;
        distB;
        xp(5) + xp(6);
        xp(5);
    ];
    R = diag([estConst.DistNoiseA,estConst.DistNoiseB, ...
        estConst.GyroNoise, estConst.CompassNoise]);
end

K = (Pp*Hk')/(Hk*Pp*Hk'+ Mk*R*Mk');
xm = xp + K*(zk - hk_xp);
Pm = (eye(6) - K*Hk)*Pp;
    
estState.xm = xm;
estState.Pm = Pm;

% Get resulting estimates and variances
% Output quantities
posEst = estState.xm(1:2);
linVelEst = estState.xm(3:4);
oriEst = estState.xm(5);
driftEst = estState.xm(6);

pxVar = estState.Pm(1,1);
pyVar = estState.Pm(2,2);
sxVar = estState.Pm(3,3);
syVar = estState.Pm(4,4);
posVar = [pxVar pyVar];
linVelVar = [sxVar syVar];
oriVar = estState.Pm(5,5);
driftVar = estState.Pm(6,6);

end

function [dxdt] = updateDiffEq(t, xpAndPp, u, estConst)
    cd = estConst.dragCoefficient;
    cr = estConst.rudderCoefficient;
    Qd = estConst.DragNoise;
    Qr = estConst.RudderNoise;
    Qb = estConst.GyroDriftNoise;
    Q = diag([Qd,Qr,Qb]);
    
    nx = 6;
    dxdt = [zeros(nx,1); zeros(nx*nx,1)];
    % q:
    dxdt(1) = xpAndPp(3); 
    dxdt(2) = xpAndPp(4);
    dxdt(3) = cos(xpAndPp(5))*(tanh(u(1))-cd*(xpAndPp(3)^2+xpAndPp(4)^2));
    dxdt(4) = sin(xpAndPp(5))*(tanh(u(1))-cd*(xpAndPp(3)^2+xpAndPp(4)^2));
    dxdt(5) = cr*u(2);
    dxdt(6) = 0;
    % matrix diff. eq:
    A_mat = A(xpAndPp,estConst,u);
    L_mat = L(xpAndPp,estConst,u);
    Pp = reshape(xpAndPp(7:42),[6,6]);
    matrix_eq = A_mat*Pp + Pp*A_mat' + L_mat*Q*L_mat';
    dxdt(7:42) = reshape(matrix_eq, [36,1]);
end

function [A_mat] = A(x,estConst,u) % (6x6)
    cd = estConst.dragCoefficient;
    cr = estConst.rudderCoefficient;
    
    A_mat = [
        0 0                    1                    0                                          0 0;
        0 0                    0                    1                                          0 0;
        0 0 -2*cd*cos(x(5))*x(3) -2*cd*cos(x(5))*x(4) -sin(x(5))*(tanh(u(1))-cd*(x(3)^2+x(2)^2)) 0;
        0 0 -2*cd*sin(x(5))*x(3) -2*cd*sin(x(5))*x(4)  cos(x(5))*(tanh(u(1))-cd*(x(3)^2+x(2)^2)) 0;
        0 0                    0                    0                                          0 0;
        0 0                    0                    0                                          0 0
    ];
end

function [L_mat] = L(x,estConst,u) % (6x3)
    cd = estConst.dragCoefficient;
    cr = estConst.rudderCoefficient;
    
    L_mat = [
                                0       0 0;
                                0       0 0;
    -cos(x(5))*cd*(x(3)^2+x(4)^2)       0 0;
    -sin(x(5))*cd*(x(3)^2+x(4)^2)       0 0;
                                0 cr*u(1) 0;
                                0       0 1
    ];
end