clc 
clear all
load heart_scale
% % 参数寻优
% [bestacc,bestc,bestg] = SVMcgForClass(heart_scale_label,heart_scale_inst,-2,4,-4,4,5,0.5,0.5,0.9);


% SVDD测试 
mu = [0,0]';
sigma = diag([1,1]);
%R=mvnrnd(mu,sigma,m) 生成200个n维正态分布数据。
data = mvnrnd(mu,sigma,200);  

trainData = data(1:100, :);
testData = data(101:200, :);


label = ones(100 ,1);
model = libsvmtrain( label , trainData , '-s 5 -t 2 -c 0.1 ' ); 
resLabel = libsvmpredict( label, testData , model );
