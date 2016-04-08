function [FT,power_produced,power_required,T_produced,T_required] = thrust(x)

    % This method of calculating thrust uses blade element theory
    %x consists of the number of blades in each rotor, rotor radius, 
    %angular velocity, number of rotors, and blade pitch, fuel consumption rate.
%     keyboard
    N = x(1);
    R = x(2)/10;
    omega = x(3)*100;
    n = x(4);
    theta = x(5)/10;
    f_rate = x(6)/10^3;

    V_c = 0;                % Vertical velocity is zero because of hover.
    a = .1091;               % Lift curve slope for naca 0012
    rho = 1.225;            % Air density, kg/m^3
    C = .09*R;              % Chord length, based on rc rotor blades (extrapolated).
    T = 0;
    Q = 0;
    P = 0;
    n_engine = 0.15;                            % Efficiency of the engine
    E_rhof = 32.4e6;                            %Energy density of fuel (Joules/L)
    f_cap = 19;             % In liters.
    W_motor = 2.1;                              %Motor weight
    W_f = f_cap*.75;                            %Fuel weight (kg)
    payload = 5;                                %Payload weight
    W_frame = 7;                                %Frame weight
    AUW = W_motor + W_f + W_frame + payload;    %All up weight in kg (about 52 lbs)
    
    m = 10;
    for i=1:m
        dr = R/m;
        r = i*R/m-.5*R/m;
        s = n*C/pi*r;
        phi = 1/16*(-sqrt(a)*sqrt(s)*sqrt(a*s+32*theta)-a*s); %atan(U_P/U_T);          % Small angle approximation because U_P << U_T.
        alpha = theta - phi;
        C_l = -98.031*alpha^3 + 0.9481*alpha^2 + 8.2078*alpha - 0.009;          % Lift coefficient based on Xfoil data.
        C_d = 51.031*alpha^4 - 0.9405*alpha^3 - 0.5316*alpha^2 + 0.0191*alpha + 0.0168;   % Drag coefficient based on xfoil data.
        t_c = C_l/6;
        lambda_i = sqrt(s*t_c/2);
        v_i = omega*r*lambda_i;
        U_T = omega*r;          % Tangential component of velocity.
        U_P = V_c + v_i;          % Perpendicular component of velocity.
  
        L = .5*rho*C*a*(U_T^2*theta-U_P*U_T);
        D = .5*rho*U_T^2*C*C_d;
%         dT = N*L*dr;
        dT = .5*rho*a*N*omega^2.*C.*(theta-phi).*r.^2.*dr;
%         dQ = N*(L*phi+D)*r*dr;
%         dQ = .5*rho*omega^2*r^3*C*(C_d+phi*C_l)*dr;
        dQ = .5*rho*omega^2*r^3*C*(C_d+phi*C_l)*dr;
        dP = dQ*omega;
        T = dT*n + T;
        Q = dQ*n + Q;
        P = dP*n + P; 
    end
    
    T_produced = T*n;
    power_required = P*n;
    power_produced = f_rate*E_rhof*n_engine;
    T_required     = AUW*9.8;   % 1.5 multiplier for controllability.
    FT = -f_cap/f_rate;             % multiply by negative to allow for minimization.
%     persistent iter;
%     if isempty(iter)
%         iter = 0;
%     end
%     iter = iter +1;
%     if iter > 1000
%         keyboard
%     end
end








