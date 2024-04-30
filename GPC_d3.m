  %2019-03-14
  
%function：基本的GPC
%parameter: 
function [y,ud_space,Theta_His] = GPC_d3(pN,PA,PB,na,nb,N,Nu,lambda,soft_ele,Cn)
%############################## 用户编辑区 ##############################
[First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,N);
%#######################################################################

%--------------------------------------------------1 表示已知量
u=zeros(1, Last_t+1-First_t);    %控制量
ud_space=zeros(1, Last_t+1-First_t);   %控制增量的存储
y=zeros(1,Last_t+1-First_t);    %输出量
w=rand(1,1+Last_t-First_t)*2-1;   %扰动。幅度范围（-1,1）的随机序列
for k = 2-First_t : 1+Last_t-First_t
    w(k) = w(k-1) + w(k);
end                             %表示ω(k)/?

%从A、B获取G
%A=[1 -0.7 0];B=[0.9 -0.4];   %第一阶段被控对象务必为此
rho=1;

%----------------------------获取多项式P、Q
M=N+Nu+(na+1)+nb;   %矩阵阶数
Theta = ones(M,1)*0.3;   %θ
Theta_His = 0.3*ones(M,1+Last_t-First_t);  %存储θ的轨迹，用于学习研究
X = zeros(M,1);    %X'=[y(t),...,y(t-(N-1)),λ△u(t+Nu-1-N),...,λ△u(t-N)
                       %-y(t-N),...,-y(t-na-N),-△u(t-N-1),...,-△u(t-N-nb)];
Ud = zeros(N+nb,1);  %[△u(t-1),...,△u(t-N-nb)]T
ud = 0;
P_theta = ones(M,M);
P_theta = P_theta+eye(M);  %随机的正定矩阵P
X_taddN = zeros(M,1);   %X(t+N), 用于求△u(t)时与θ做矩阵乘法

tic;
%31-First_t
for k = 1-First_t : 1+Last_t-First_t   %时刻指针k, k=1-First_t对应实际时刻t=0
%for t = 0         :  Last_t
    %----------------------------------------3.1 检测k时刻输出y(k)
    y(k) = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn);    %获取y(k)
    %if(k==1-First_t+50 || k==1-First_t+130)  %添加人工干扰
        %y(k)=y(k)+10;
    %end
    %----------------------------------------3.2 参数估计
    %------------------------------3.2.1 更新Ud,X
    Ud = [Ud(1:N+nb-1,1)];
    Ud = [ud;Ud];
    X = [X(1:M-1,1)];
    X = [y(k);X];
    X(N+1) = lambda*Ud(N+1-Nu);
    X(N+Nu+1) = -y(k-N);
    X(N+Nu+na+2) = -Ud(N+1);
    %-----------------------------3.2.2 更新θ
    epsilon = Ud(N)-X'*Theta;
    Mid = P_theta*X;         %P(t-1)*X(t)
    Mid = Mid/(X'*Mid+rho);  %P(t-1)*X(t)/(rho+X(t)'*P(t-1)*X(t))
    Theta = Theta+Mid*epsilon;  
    P_theta = (P_theta-Mid*X'*P_theta)/rho;
    Theta_His(:,k) = Theta;  %记录θ轨迹
    %----------------------------------------3.3 实时更新控制率
    %求出△u(t),u(t)
    for i=M:-1:N+Nu+2
        X_taddN(i) = X_taddN(i-1);
    end
    X_taddN(N+Nu+1) = -y(k);  %历史输出项
    X_taddN(N+Nu+na+2) = -Ud(1);  %历史输入项
    %----- 设定值项
    mid = yr(k)*(1-soft_ele);
    X_taddN(N)=y(k)*soft_ele + mid;
    for i=N-1:-1:1
        X_taddN(i) = X_taddN(i+1)*soft_ele + mid;
    end
    
    ud = X_taddN'*Theta;
    u(k) = u(k-1) + ud;
    if u(k)<u_min           %给u(t)一些限制条件
        u(k) = u_min;
        ud = u_min - u(k-1);
    elseif(u(k)>u_max)
        u(k) = u_max;
        ud = u_max - u(k-1);
    end
    ud_space(k)=ud;
end
run_time=toc;

figure(5);
uicontrol('Style','text','Position',[10 0 150 20],'String',['运行时间: ',num2str(run_time),'s']);
        %运行时间显示

subplot(211);
plot(t,yr,'b--',t,y,'r-');
grid on;
title('简单的 广义预测自适应控制 直接算法仿真');
%xlim([First_t Last_t]);
axis([First_t,Last_t,-10,60]);
ylabel('输出y(k)');
legend('设定值yr','实际值y','Location','Best');

subplot(212);
plot(t,ud_space);
grid on;
xlim([First_t Last_t]);
%axis([First_t,Last_t,-40,40]);
ylabel('控制增量u_d(k)');
xlabel('k')

%--------------------------------------------------5 编辑返回值
y = [y(1,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
u = [u(1,1-First_t:1+Last_t-First_t)];
ud_space = [ud_space(1,1-First_t:1+Last_t-First_t)];
Theta_His = [Theta_His(:,1-First_t:1+Last_t-First_t)];
DataTXT = fopen('.\GPC_simulate_data\GPC_d3__data.txt','w');%打开txt文件
fprintf(DataTXT,'%s\t%s\t%s\t\t%s\r\n','时刻k','输出y','控制增量ud','各项估计参数θ(1~?)');
data_len = length(y(:));
theta_len = length(Theta_His(:,1));
for j=1:data_len
    fprintf(DataTXT,'%d\t',j-1);
    fprintf(DataTXT,'%.3f\t%.3f\t\t',y(j),ud_space(j));
    for i = 1:theta_len
        fprintf(DataTXT,'%.4f\t',Theta_His(i,j));
    end
    fprintf(DataTXT,'\r\n');
end
fclose(DataTXT);
