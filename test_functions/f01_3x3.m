f  =@(x) x(1).^4+3*x(1).^3-2*x(1).^2+x(2).^4+9*x(2).^2+x(1).*x(2)+x(3).^2;
df =@(x)[4*x(1).^3+9*x(1).^2-4*x(1)+x(2);4*x(2).^3+18*x(2)+x(1);2*x(3)];
hf =@(x)[12*x(1).^2+18*x(1)-4,1,0;1,12*x(2).^2+18,0;0,0,2];
x0 = [-3;-5;1];