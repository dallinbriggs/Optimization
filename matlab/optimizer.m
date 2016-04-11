function [opt_A,fval,exitflag,output_fmin,lambda,gradient,constraints] = optimizer()



options = optimoptions('fmincon',...
        'Algorithm', 'sqp', ...  % choose one of: 'interior-point', 'sqp', 'active-set', 'trust-region-reflective'
        'AlwaysHonorConstraints', 'bounds', ...  % forces interior point algorithm to always honor bounds
        'display', 'iter-detailed', ...  % display more information
        'Tolx', 1e-15, ...
        'MaxIter', 10000, ...  % maximum number of iterations
        'MaxFunEvals', 100000, ...  % maximum number of function calls
        'TolCon', 1e-7, ...  % convergence tolerance on constraints
        'TolFun', 1e-7, ...  % convergence tolerance on function value
        'GradObj', 'on', ...  % supply gradients of objective
        'GradConstr', 'off', ...  % supply gradients of constraints
        'Diagnostics', 'on', ...  % display diagnotic information
        'PlotFcns',{@optimplotfval,@optimplotfirstorderopt});

A0 = [2,.5*10,2*pi*2000/60/100,4,5*pi/180*10,.001*10^3]; % Initial guess on number of blades in each rotor, 
                                       %rotor radius, angular velocity, number of rotors, 
                                       %and blade pitch, fuel consumption rate.                

Aeq = [];
Beq = [];
lb = [2, .01*10, 2*pi*100/60/100, 1, 0, 0.0001848346*10^3]; 
ub = [2, .7*10, 2*pi*4500/60/100, 4, 20*pi/180*10, 0.0246446*10^3];


[opt_A,fval,exitflag,output_fmin,lambda,gradient] = fmincon(@obj, A0, [], [], Aeq, Beq, lb, ub,@con,options);

[~,constraints.power_produced, constraints.power_required, constraints.thrust_produced, ...
    constraints.thrust_required] = thrust(opt_A);


%     ---------- Objective Function ------------------
    function [J,g] = obj(x)
        
        [f,~,~,~,~] = thrust(x);
        dfdx = grad(x);
        J = f / 1e4;  % scale to make objective of order(1)
        g = dfdx / 1;
    end



    function [cin, ceq, dcin, dceq] = con(x)
        
        [~,power_produced,power_required,thrust_produced,thrust_required] = thrust(x);
   
        cin(1) = 1*(power_required - power_produced);  % of the form c(x) <= 0
        cin(2) = 1*(thrust_required - thrust_produced);
        ceq = [];
        dcin = [];  
        dceq = [];
 
    end

end

