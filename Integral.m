% It is a code to calculate the integral
% "var" is a variable which that be integrated
% "var" should be a vector
% "dt" is a time interval between each variable
% "dt" should be a constant
% "int_v" is an output vector which are integrated values of "var"
% Copyright (c) 2023 Hong Gi Yeom

function [int_v]=Integral(var, dt)

int_v=zeros(1,length(var));
int_v(1)=var(1)*dt;
for k=2:1:length(var);
    int_v(k)=int_v(k-1)+(var(k)*dt);
end
