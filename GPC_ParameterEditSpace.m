%2019-04-16

%function：生成可供用户编辑的参数，用户可在此处对软件进行更高级的参数配置
%parameter: 参数含义与GPC各算法中的含义均相同
function [First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,N)
%############################## 用户编辑区 ##############################
%---------- 时域范围
First_t = -10;
if (abs(First_t)<=N)
    First_t = -N - 5;
end
Last_t=-1;
tim_len = 80;  %每个被控对象所占的时刻数
for i=1:pN
    Last_t = Last_t + tim_len; %79,159,...
end
t = First_t : Last_t;  %-10 ~ 79/159/239
%===================
%---------- 设定值/预测值形式
yr=zeros(1, Last_t+1-First_t);
yr = yr+10*stepfun(t,0);
for i=1:pN
    %yr = yr + 1*stepfun(t,(i-1)*tim_len+30);
    %yr = yr - 1*stepfun(t,(i-1)*tim_len+60);
    yr = yr + 10*stepfun(t,(i-1)*tim_len+20);%段内延迟20的位置，有幅度为10的阶跃点
    yr = yr + 10*stepfun(t,(i-1)*tim_len+40);%段内延迟40，阶跃幅度10
    yr = yr - 20*stepfun(t,(i-1)*tim_len+60);%段内延迟60，阶跃幅度-20
end  %[0,79] [80 159] [160 239]
%===========================
u_max=999999999999999;  %控制量上限，即控制量的约束条件
u_min=-999999999999999;   %控制量下限
%#######################################################################
