  %2019-03-14
  
%function：GPC自适应算法
%parameter: 
%return: Theta_His: 估计参数的轨迹
%note:自适应算法，ρ取小一点效果好，过小稳定性差
%     A、B初值对系统无影响
function [y,ud_space,Theta_His] = GPC_2(pN,PA,PB,na,nb,N1,Nu,lambda,rho,soft_ele,Cn)
%############################## 用户编辑区 ##############################
[First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,N1);
%#######################################################################

%--------------------------------------------------2 表示其他已知量
u=zeros(1, Last_t+1-First_t);    %控制量
ud_space=zeros(1, Last_t+1-First_t);   %控制增量的存储
y=zeros(1,Last_t+1-First_t);    %输出量
w=rand(1,1+Last_t-First_t)*2-1;   %扰动。幅度范围（-1,1）的随机序列
for k = 2-First_t : 1+Last_t-First_t
    w(k) = w(k-1) + w(k);
end                             %表示ω(k)/?

A = ones(1,na+1);   %A(1)永远等于1
B = ones(1,nb+1);
n_theta = na+nb+1;                 %产生θ的阶数
Theta = [-3;0;0.5;0];        %θ的初值
  %θ各项表示：
Theta_His = ones(n_theta,1+Last_t-First_t)*1;
P_theta = ones(n_theta,n_theta)*2;
P_theta = P_theta+eye(n_theta)*1;  %随机的正定矩阵P
X = zeros(n_theta,1);          %微分项向量[1~na,na+1~n_theta]

ud = 0;                        %用于存储△u(t)
yd = 0;                        %用于存储△y(t)
Ud = zeros(1,nb); Ud_p=1;  %存储[△u(t-1),...,△u(t-nb)]
Yr = zeros(N1,1);          %存储柔化后的设定值项

%--------------------------------------------------3 控制过程
tic;
for k = 1-First_t : 1+Last_t-First_t   %时刻指针t, t_p=1-First_t对应实际时刻t=0
%for t = 0         :  Last_t
    %------------------------------3.1 检测t_p时刻输出
    y(k) = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn);    %获取y(k)
    %if(k==1-First_t+50 || k==1-First_t+130)  %添加人工干扰
    %    y(k)=y(k)+10;
    %end
    %------------------------------3.2 参数估计
    for i=n_theta:-1:2
         X(i) = X(i-1);  %更新X
    end
    X(na+1) = ud;  %△u(t-1)进入X
    X(1) = -yd;    %-△y(t-1)进入X
    yd = y(k)-y(k-1); %update △y(t)
    %-------------------- update θ(t),P(t-1)
    Mid = P_theta*X;         %P(t-1)*X(t)
    Mid = Mid/(X'*Mid+rho);  %P(t-1)*X(t)/(rho+X(t)'*P(t-1)*X(t))
    Theta = Theta+Mid*(yd-X'*Theta);  
    P_theta = (P_theta-Mid*X'*P_theta)/rho;
    Theta_His(:,k) = Theta;
    %------------------------------3.3 实时更新控制率
    for i=1:na
        A(i+1) = Theta(i);
    end
    for i=1:nb+1
        B(i) = Theta(i+na);
    end
    [P,Alpha,Beta]=GPC_getCtrlRule(A,na,B,nb,N1,Nu,lambda);
    %------------------------------3.4 更新控制量△u(t)
    ud = 0;    %△u(t)
    
    mid = yr(k)*(1-soft_ele);
    Yr(1) = y(k)*soft_ele + mid;
    for j=2:N1
        Yr(j) = Yr(j-1)*soft_ele +mid;
    end
    ud = ud + P*Yr;    %其中P是横向量
    
    for j=1:na+1
        ud = ud - Alpha(j)*y(k+1-j);
    end
    
    for j=1:nb
        ud = ud - Beta(j)*Ud(mod(Ud_p-2+j,nb)+1);  %Ud_p指向△u(t-1),Ud_p+1指向△u(t-2)
    end    %已求出△u(t)
    
    %----------2 求u(k)
    u(k) = u(k-1) + ud;
    if u(k)<u_min           %给u(t)一些限制条件
        u(k) = u_min;
        ud = u_min - u(k-1);
    elseif(u(k)>u_max)
        u(k) = u_max;
        ud = u_max - u(k-1);
    end
    ud_space(k)=ud;
    %---------- 更新△u(k-1)存储数组
    Ud_p = Ud_p-1;
    if Ud_p < 1
        Ud_p = nb;
    end
    Ud(Ud_p) = ud;     %以上5行。更新△u(t-1)
end
run_time=toc;

%--------------------------------------------------4 绘图
figure(4);
uicontrol('Style','text','Position',[10 0 150 20],'String',['运行时间: ',num2str(run_time),'s']);
        %运行时间显示
subplot(211);
plot(t,yr,'b--',t,y,'r-');
legend('设定值yr','实际值y','Location','Best');
axis([First_t Last_t -10 50]);
%xlim([First_t Last_t]);   %只限制横坐标范围
grid on;
title('广义预测自适应控制算法仿真');
ylabel('输出y(k)');
subplot(212);
plot(t,ud_space);
%plot(t,Theta_His(1,:));
%axis([First_t Last_t -10 100]);
xlim([First_t Last_t]);
grid on;
ylabel('控制增量u_d(k)');
%ylabel('a_1(k)');
xlabel('k')

%--------------------------------------------------5 编辑返回值
%------------------------------ 整理历史轨迹
y = [y(1,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
u = [u(1,1-First_t:1+Last_t-First_t)];
ud_space = [ud_space(1,1-First_t:1+Last_t-First_t)];
Theta_His = [Theta_His(:,1-First_t:1+Last_t-First_t)];

%------------------------------ 历史轨迹输出至txt
DataTXT = fopen('.\GPC_simulate_data\GPC_2__data.txt','w'); %打开txt文件
fprintf(DataTXT,'%s\t%s\t%s\t\t%s\r\n','时刻k','输出y','控制增量ud','各项估计参数θ(1~?)');
    %输入标题
data_len = length(y(:));  %数据长度（仿真时间）获取
theta_len = length(Theta_His(:,1));  %参数个数获取
for j=1:data_len
    fprintf(DataTXT,'%d\t',j-1);
    fprintf(DataTXT,'%.3f\t%.3f\t\t',y(j),ud_space(j));
    for i = 1:theta_len
        fprintf(DataTXT,'%.4f\t',Theta_His(i,j));
    end
    fprintf(DataTXT,'\r\n');
end

%------------------------------ 求取对应PID参数（最后时刻）
[E,F,G,H] = GPC_getEFGH(A,na,B,nb,N1,Nu);
Kp = 0;
Ki = 0;
Kd = 0;
for j=1:N1
    Kp = Kp - P(j)*(F(j,2)+2*F(j,3));
    Ki = Ki + P(j);
    Kd = Kd + P(j)*F(j,3);
end
fprintf('Kp = %.4f\n',Kp)
fprintf('Ki = %.4f\n',Ki)
fprintf('Kd = %.4f\n\n',Kd)

fclose(DataTXT);
