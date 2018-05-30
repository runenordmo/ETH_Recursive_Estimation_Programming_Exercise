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

% PRIOR UPDATE

% Set up differential equation for prior update Hybrid EKF
nx = 6;
xpAndPp0 = [estState.xm; reshape(estState.Pm,[nx*nx 1])];
[t,xpAndPp] = ode45(@(t,xpAndPp) priorUpdateDiffEq(t,xpAndPp,actuate,estConst),[tm-dt tm],[xpAndPp0]);
xpAndPpEnd = xpAndPp(end,:)';

% Extract prior estimate and variance
xp = xpAndPpEnd(1:nx);
Pp = reshape(xpAndPpEnd(nx+1:nx+nx*nx),[nx nx]);

% MEASUREMENT UPDATE

% Set posterior equal to prior initially
xm = xp; Pm = Pp;

% Extract measurement
zk = sense';

% Perform measurement update according to availability of measurement 3 
isMeasurement3Available = isfinite(zk(3));
if(isMeasurement3Available)
    zk_available = zk;
else
    zk_available = [zk(1:2);zk(4:5)];
end
[H,M,h,R] = getMeasurementModel(xp,estConst,isMeasurement3Available);

% Update Kalman gain and posterior mean and variance
K = (Pp*H')/(H*Pp*H'+ M*R*M');
xm = xp + K*(zk_available - h);
Pm = (eye(6) - K*H)*Pp;

% Set posterior mean and variance in return struct estState
estState.xm = xm;
estState.Pm = Pm;

% Get resulting estimates and variances
% Resulting estimates
posEst = xm(1:2);
linVelEst = xm(3:4);
oriEst = xm(5);
driftEst = xm(6);

% Resulting variances
posVar = [Pm(1,1) Pm(2,2)];
linVelVar = [Pm(3,3) Pm(4,4)];
oriVar = Pm(5,5);
driftVar = Pm(6,6);

end

function [dxpAndPp_dt] = priorUpdateDiffEq(t, xpAndPp, u, estConst)
    nx = 6;
    dxpAndPp_dt = [zeros(nx,1); zeros(nx*nx,1)];
    
    % Process equation for the mean:
    dxpAndPp_dt(1:nx) = q(xpAndPp, estConst,u);
    % Matrix differential equation for variance:
    Pp = reshape(xpAndPp(nx+1:nx+nx*nx),[nx,nx]);
    A_mat = A(xpAndPp,estConst,u);
    L_mat = L(xpAndPp,estConst,u);
    Q = diag([estConst.DragNoise,estConst.RudderNoise,estConst.GyroDriftNoise]);
    
    matrix_eq = A_mat*Pp + Pp*A_mat' + L_mat*Q*L_mat';
    dxpAndPp_dt(nx+1:nx+nx*nx) = reshape(matrix_eq, [nx*nx,1]);
end

% Process dynamics
function [q_vec] = q(x,estConst,u) % (6x1)
    cd = estConst.dragCoefficient;
    cr = estConst.rudderCoefficient;
    
    q_vec = [
        x(3); 
        x(4);
        cos(x(5))*(tanh(u(1))-cd*(x(3)^2+x(4)^2));
        sin(x(5))*(tanh(u(1))-cd*(x(3)^2+x(4)^2));
        cr*u(2);
        0;
    ];
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

function [H,M,h,R] = getMeasurementModel(xp,estConst,isMeasurement3Available)
    % Set up convenient variables
    xDistA = xp(1) - estConst.pos_radioA(1);
    yDistA = xp(2) - estConst.pos_radioA(2);
    xDistB = xp(1) - estConst.pos_radioB(1);
    yDistB = xp(2) - estConst.pos_radioB(2);
    xDistC = xp(1) - estConst.pos_radioC(1);
    yDistC = xp(2) - estConst.pos_radioC(2);
    distA = sqrt(xDistA^2 + yDistA^2);
    distB = sqrt(xDistB^2 + yDistB^2);
    distC = sqrt(xDistC^2 + yDistC^2);

    % Get actual measurement model
    if(isMeasurement3Available)
        H = [ ...
            xDistA/distA yDistA/distA 0 0 0 0;
            xDistB/distB yDistB/distB 0 0 0 0;
            xDistC/distC yDistC/distC 0 0 0 0;
                       0            0 0 0 1 1;
                       0            0 0 0 1 0;
        ];
        M = eye(5);
        h = [ ...
            distA;
            distB;
            distC;
            xp(5) + xp(6);
            xp(5);
        ];
        R = diag([estConst.DistNoiseA,estConst.DistNoiseB, ...
            estConst.DistNoiseC, estConst.GyroNoise, estConst.CompassNoise]);
    else
        H = [ ...
            xDistA/distA yDistA/distA 0 0 0 0;
            xDistB/distB yDistB/distB 0 0 0 0;
                       0            0 0 0 1 1;
                       0            0 0 0 1 0;
        ];
        M = eye(4);
        h = [ ...
            distA;
            distB;
            xp(5) + xp(6);
            xp(5);
        ];
        R = diag([estConst.DistNoiseA,estConst.DistNoiseB, ...
            estConst.GyroNoise, estConst.CompassNoise]);
    end
    
end