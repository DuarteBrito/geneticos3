%y = y_basis + sum(ai * fi)

x = 0:0.1:1;
t1s = linspace(0,1,11);
%t1s = [0.33 0.66];
%t1 = 0.66;
t2 = 5;

y = x;

%fi
for t1 = t1s
y = x;
%for ii = 1:length(y)
    %y(ii) = sin(pi*x(ii)^(log(0.5)/log(t1)))^t2;
    
%y = f(x);
%end
y = sin(pi*x.^(log(0.5)/log(t1))).^t2;
plot(x,y)
hold on
end

% hold on
% plot(x,y)

