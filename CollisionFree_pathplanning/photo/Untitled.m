clc
close all
clear all
x=1:10;%give your x here
h=1:7; %give your h here
m=length(x);
n=length(h);
X=[ones(1,n)];
H=[h,zeros(1,m)];
%convolution


Y = conv(h,h);


stem(Y);
ylabel('Y[n]');
xlabel('n');