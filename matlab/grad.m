function [g] = grad(x0)

addpath(strcat(pwd,'/autodiff'));

x = adiff(x0);
[a,~,~,~,~] = thrust(x);
[~,g] = adiffget(a);


end