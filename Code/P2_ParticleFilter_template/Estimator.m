function [postParticles] = Estimator(prevPostParticles, sens, act, init)
% [postParticles] = Estimator(prevPostParticles, sens, act, init)
%
% The estimator function. The function will be called in two different
% modes: If init==1, the estimator is initialized. If init == 0, the
% estimator does an iteration for a single sample time interval Ts (KC.ts)
% using the previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurements and control inputs.
%
% You must edit this function.
%
% Inputs:
%   prevPostParticles   previous posterior particles at discrete time k-1,
%                       which corresponds to continuous time t = (k-1)*Ts
%                       The variable is a struct whose fields are arrays
%                       that correspond to the posterior particle states.
%                       The fields are: (N is number of particles)
%                       .x = (2xN) array with the x-locations (metres)
%                       .y = (2xN) array with the y-locations (metres)
%                       .h = (2xN) array with the headings (radians)
%                       The first row in the arrays corresponds to robot A.
%                       The second row corresponds to robot B.
%
%   sens                Sensor measurements at discrete time k (t = k*Ts),
%                       [4x1]-array, an Inf entry indicates no measurement
%                       of the corresponding sensor.
%                       sens(1): distance reported by sensor 1 (metres)
%                       sens(2): distance reported by sensor 2 (metres)
%                       sens(3): distance reported by sensor 3 (metres)
%                       sens(4): distance reported by sensor 4 (metres)
%
%   act                 Control inputs u at discrete time k-1, which are
%                       constant during a time interval Ts:
%                       u(t) = u(k-1) for (k-1)*Ts <= t < k*Ts
%                       [2x1]-array:
%                       act(1): velocity of robot A, u_A(k-1) (metres/second)
%                       act(2): velocity of robot B, u_B(k-1) (metres/second)
%
%   init                Boolean variable indicating wheter the estimator
%                       should be initialized (init = 1) or if a regular
%                       estimator update should be performed (init = 0).
%                       OPTIONAL ARGUMENT. By default, init = 0.
%
% Outputs:
%   postParticles       Posterior particles at discrete time k, which
%                       corresponds to the continuous time t = k*Ts.
%                       The variable is a struct whose fields are arrays
%                       that correspond to the posterior particle states.
%                       The fields are: (N is number of particles)
%                       .x = (2xN) array with the x-locations (metres)
%                       .y = (2xN) array with the y-locations (metres)
%                       .h = (2xN) array with the headings (radians)
%                       The first row in the arrays corresponds to robot A.
%                       The second row corresponds to robot B.
%
% Class:
% Recursive Estimation
% Spring 2018
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control

% Check if init argument was passed to estimator:
if(nargin < 4)
    % if not, set to default value:
    init = 0;
end

%% Mode 1: Initialization
% Set number of particles:
%global Nparticles;
%N = Nparticles;
N = 10000; % obviously, you will need more particles than 10.
if (init)
    % Do the initialization of your estimator here!
    % These particles are the posterior particles at discrete time k = 0
    % which will be fed into your estimator again at k = 1
    % Replace the following:
    postParticles.x = zeros(2,N);
    postParticles.y = zeros(2,N);
    postParticles.h = zeros(2,N);
    
    % Initialize A
    postParticles.x(1,:) = 2 * KC.L;
    for i = 1:1:N
        randa = rand(1);
        if randa <= 0.5
            postParticles.h(1,i) = pi/2 + rand(1) * pi/2;
        else
            postParticles.y(1,i) = KC.L;
            postParticles.h(1,i) = -pi/2 - rand(1) * pi/2;
        end
    end
    
    % Initialize B
    for i = 1:1:N
        randb = rand(1);
        if randb <= 0.5
            postParticles.h(2,i) = rand(1) * pi/2;
        else
            postParticles.y(2,i) = KC.L;
            postParticles.h(2,i) = rand(1) * -pi/2;
        end
    end
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If init = 0, we perform a regular update of the estimator.

% Implement your estimator here!

% Replace the following:
postParticles.x = zeros(2,N);
postParticles.y = zeros(2,N);
postParticles.h = zeros(2,N);

% Step 1 : Prior Update
xp = zeros(2,N);
yp = zeros(2,N);
hp = zeros(2,N);
postParticles.h = prevPostParticles.h;

for i = 1:1:N
    vs_a = tri(KC.vsbar,1);
    ua = act(1) * (1 + vs_a);
    xp_a = ua * cos(prevPostParticles.h(1,i)) * KC.ts + prevPostParticles.x(1,i);
    yp_a = ua * sin(prevPostParticles.h(1,i)) * KC.ts + prevPostParticles.y(1,i);
    [xp_a, yp_a, hp_a] = check_bounce(prevPostParticles.x(1,i),...
        prevPostParticles.y(1,i), xp_a, yp_a, ua, KC.ts, prevPostParticles.h(1,i));
    xp(1,i) = xp_a;
    yp(1,i) = yp_a;
    hp(1,i) = hp_a;
    
    vs_b = tri(KC.vsbar,1);
    ub = act(2) * (1 + vs_b);
    xp_b = ub * cos(prevPostParticles.h(2,i)) * KC.ts + prevPostParticles.x(2,i);
    yp_b = ub * sin(prevPostParticles.h(2,i)) * KC.ts + prevPostParticles.y(2,i);
    [xp_b, yp_b, hp_b] = check_bounce(prevPostParticles.x(2,i),...
        prevPostParticles.y(2,i), xp_b, yp_b, ub, KC.ts, prevPostParticles.h(2,i));
    xp(2,i) = xp_b;
    yp(2,i) = yp_b;
    hp(2,i) = hp_b;
end
hp = mod(hp+pi,2*pi)-pi;
%% ------------------------------------------------------------------------
% Step 2 : Measurement Update

sensPos = [2*KC.L 2*KC.L  0      0;
           0      KC.L    KC.L   0];
anyMeasurment = false;
Betas = ones(1,N);
for i = 1:4
    sensorA = i < 3;
    if (isfinite(sens(i)))
        anyMeasurment = true;
        % Apriori estimate of distance
        d_A = vecnorm([xp(1,:);yp(1,:)] - sensPos(:,i));
        d_B = vecnorm([xp(2,:);yp(2,:)] - sensPos(:,i));
        % likelihood of measurements. Using an uniform triangular
        % assumption of measurement noise
        wbar = KC.wbar;
        if wbar < 0.05
            wbar = 0.05;
        end
        p_A = max(0,(1/wbar)*(1-abs(d_A-sens(i))/wbar));
        p_B = max(0,(1/wbar)*(1-abs(d_B-sens(i))/wbar));
        
        if sensorA
            % sbar is the prob of sensors measuring distance
            % to wrong robot.
            betas = [KC.sbar (1-KC.sbar)]*[p_B; p_A];
        else % sensorB
            betas = [KC.sbar (1-KC.sbar)]*[p_A; p_B];
        end
        
        % If all particles have zero probability, skip and
        % assume roughening will fix it.
        %if sum(betas) <= 0
        %    continue;
        %end
        Betas = Betas.*betas;
    end
end

E_x = 2*KC.L;
E_y = KC.L;
E_h = pi;

if anyMeasurment
    alpha = 1/sum(Betas);
    Betas = alpha*Betas;
    Beta_cum = cumsum(Betas);
    numericalError = any(isnan(Beta_cum));
    if not(numericalError)
        % Resample
        for n = 1:N
           randNumber = rand;
           ind = find(Beta_cum >= randNumber,1,'first');
           postParticles.x(:,n) = xp(:,ind);
           postParticles.y(:,n) = yp(:,ind);
           postParticles.h(:,n) = hp(:,ind);
        end
        % Use fine Roughening
        E_x = max(postParticles.x,[],2)-min(postParticles.x,[],2);
        E_y = max(postParticles.y,[],2)-min(postParticles.y,[],2);
        E_h = max(postParticles.h,[],2)-min(postParticles.h,[],2);
    else % Numerical error (all particles has zero probability)
         % Use prior estimate
        postParticles.x = xp;
        postParticles.y = yp;
        postParticles.h = hp;
        % Use rough Roughening
    end
else % No new measurments
    % Use prior estimate
    postParticles.x = xp;
    postParticles.y = yp;
    postParticles.h = hp;
    % Use fine Roughening
    E_x = max(postParticles.x,[],2)-min(postParticles.x,[],2);
    E_y = max(postParticles.y,[],2)-min(postParticles.y,[],2);
    E_h = max(postParticles.h,[],2)-min(postParticles.h,[],2);
end

% Roughening
K = 0.005;

std_x = K*E_x*N^(1/6);
std_y = K*E_y*N^(1/6);
std_h = K*E_h*N^(1/6);

delta_x = transpose(std_x.^2)*randn(2,N);
delta_y = transpose(std_y.^2)*randn(2,N);
delta_h = transpose(std_h.^2)*randn(2,N);

postParticles.x = max(0, min(2*KC.L, postParticles.x + delta_x));
postParticles.y = max(0, min(KC.L, postParticles.y + delta_y));
postParticles.h = mod(postParticles.h + delta_h + pi,2*pi) - pi;

end % end estimator

%--------------------------------------------------------------------------
% Probability distrebution generators
function X = tri(b,n)
    X = zeros(n,1);
    for i = 1:n
        z = rand;
        if sqrt(2 * z * b^2) - b < 0
            X(i) = sqrt(2 * z * b^2) - b;
        else
            X(i) = b - sqrt(2 * (1-z) * b^2);
        end
    end
end

function X = quad(v,n)
    X = zeros(n,1);
    for i = 1:1:n
        sign = rand;
        if sign > 0.5
            X(i) = random('beta', 3, 1) * v;
        else
            X(i) = -random('beta', 3, 1) * v;
        end
    end
end
%--------------------------------------------------------------------------
function [x,y,h] = check_bounce(x_prev,y_prev,x,y,u,ts,h_prev)
    h = h_prev;
    if(x >= 0 && x <= 2*KC.L && y >= 0 && y <= KC.L)
        return;
    end
    
    if y < 0
        t_inter = (0 - y_prev) / (u * sin(h_prev));
        x_bounces = u * cos(h_prev)  * t_inter + x_prev;
        if(x_bounces >= 0 && x_bounces <= 2*KC.L)
            if(h_prev < -pi/2)
                h = pi - (pi + h_prev) * (1 + quad(KC.vbar,1));
            else
                h = -h_prev * (1 + quad(KC.vbar,1));
            end
            t_left =  ts- t_inter;
            x = u * cos(h) * t_left + x_bounces;
            y = u * sin(h) * t_left;
            [x, y, h] = check_bounce(x_bounces, 0, x, y, u, t_left, h);
            return;
        end
    elseif(y > KC.L)  %bounce the upper, horizontal wall
        t_inter = (KC.L - y_prev) / (u * sin(h_prev));
        x_bounces = u * cos(h_prev) * t_inter + x_prev;
        if(x_bounces >= 0 && x_bounces <= 2*KC.L)
            if(h_prev < pi/2)
                h = -h_prev * (1 + quad(KC.vbar,1));
            else
                h = -pi +  (pi - h_prev) * (1 + quad(KC.vbar,1));
            end
            t_left =  ts- t_inter;
            x = u * cos(h) * t_left + x_bounces;
            y = u * sin(h) * t_left + KC.L;
            [x, y, h] = check_bounce(x_bounces, KC.L, x, y, u, t_left, h);
            return;
        end
    end
    
    if(x < 0)     %bounce the left, vertical wall
        t_inter = (0 - x_prev) / (u * cos(h_prev));
        y_bounces = u * sin(h_prev) * t_inter + y_prev;
        if(y_bounces >= 0  && y_bounces <= KC.L)
            if(h_prev < 0)
                h = - (pi/2 +(pi/2 + h_prev) * (1 + quad(KC.vbar,1)));
            else
                h = pi/2 -(h_prev - pi/2) * (1 + quad(KC.vbar,1));
            end
            t_left =  ts- t_inter;
            x = u * cos(h) * t_left;
            y = u * sin(h) * t_left + y_bounces;
            [x, y, h] = check_bounce(0,y_bounces, x, y, u, t_left, h);
            return;
        end
    elseif(x > 2*KC.L)    %bounce the right, vertical wall
        t_inter = (2*KC.L - x_prev) / (u * cos(h_prev));
        y_bounces = u * sin(h_prev) * t_inter + y_prev;
        if(y_bounces >= 0  && y_bounces <= KC.L)
            if(h_prev < 0)
                h = -pi/2 - (h_prev + pi/2) * (1 + quad(KC.vbar,1));
            else
                h = pi/2 + (pi/2 - h_prev) * (1 + quad(KC.vbar,1));
            end
            t_left =  ts- t_inter;
            x = u * cos(h) * t_left + 2 * KC.L;
            y = u * sin(h) * t_left + y_bounces;
            [x, y, h] = check_bounce(2 * KC.L,y_bounces, x, y, u, t_left, h);
            return;
        end
    end
end
    
