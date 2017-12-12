function yi = interp1(x,y,xi,method)
% interp1(x,y,xi) does interpolation for autodiff objects 
%   - note that x must be a double

if nargin==3, method='linear'; end
switch lower(method)
case {'linear','*linear'}
   x  = x(:);
   [x,xidx] = sort(x); % which means that x can't be an adiff object
   y  = y(xidx);
   dy = diff(y)./diff(x);
   xidx = zeros(size(xi));
   for i=1:length(xi)
      if xi(i)==x(end)
         xidx(i) = length(x)-1;
      else
         xidx(i) = sum(x<=xi(i));
      end
   end
   yi = y(xidx)+dy(xidx).*(xi(:)-x(xidx));
   
otherwise
   error([method,' method is not supported for autodiff objects']);
end
