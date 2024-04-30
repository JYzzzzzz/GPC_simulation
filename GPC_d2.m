%2019-04-05

%function：改进的GPC直接算法，基于对象阶跃响应前N时刻估计值的GPC直接算法
%parameter: lambda：测试时取1.25
%           forget_ele：遗忘因子，测试时取0.95
%           soft_ele：柔化因子，测试时取0
function [y,ud_space,Theta_His] = GPC_d2(pN,PA,PB,na,nb,N,Nu,lambda,forget_ele,soft_ele,Cn)
%############################## 用户编辑区 ##############################
[First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,N);
M=0.01;    %可接受误差上限
C=2;       %正则化防出错数
%#######################################################################

%--------------------------------------------------1 表示其他已知量
u=zeros(1, Last_t+1-First_t);    %控制量
ud_space=zeros(1, Last_t+1-First_t);   %控制增量的存储
y=zeros(1,Last_t+1-First_t);    %输出量
w=rand(1,1+Last_t-First_t)*2-1;   %扰动。幅度范围（-1,1）的随机序列
for k = 2-First_t : 1+Last_t-First_t
    w(k) = w(k-1) + w(k);
end                             %表示ω(k)/?
%------------------------------1.2 计算参数P、Q
close_ele=0.5;     %用于拟定阶跃响应前N项g0~gN-1的一个参数，
                   %实际应用中，g最好使用测量值
                   %当然经过验证得，g适当的随便取值，系统也可以良好的工作
Itmd_val=close_ele*ones(N,1);   %存g0~gN-1
if(N>1)
    for i=2:N
        Itmd_val(i) = (1-Itmd_val(i-1))*close_ele+Itmd_val(i-1);
    end
end 
%Itmd_val=[0.40;1.06;1.45;1.48;1.30;1.13];
G=zeros(N,N);    %矩阵G
for row=1:N
    i=0;
    while row+i<=N
        G(row+i,1+i)=Itmd_val(row);
        i=i+1;
    end
end
G = [G(:,1:Nu)];
Q = G'*G + lambda*eye(Nu);
Q = inv(Q);
P = Q * G';
P = [P(1,:)];
Q = [Q(1,:)];

%-------------------------------------------------2 空间开辟
Yr = zeros(N,1);       %含有柔化因子影响的yr序列
Yd = zeros(na+N,1);    %存储[△y(k),...,△y(k-(na+N-1))]
Ud = zeros(nb+N,1);    %存储[△u(k-1),...,△u(k-N-nb)]
ud = 0;                %用于存储△u(k)
Theta = ones(na+nb,1)*0;    %θ的初值
Theta_His = ones(na+nb,1+Last_t-First_t)*0;
Y = zeros(na+N,1);  %正则化过程中需用到，存储[y(k),...,y(k-(na+N-1))]各项的绝对值乘以2
U = zeros(nb+N,1);  %存储[u(k-1),...,u(k-nb-N)]各项的绝对值乘以2

%--------------------------------------------------3 控制过程
tic;
for k = 1-First_t : 1+Last_t-First_t   %时刻指针t, t_p=1-First_t对应实际时刻t=0
    %------------------------------3.1 检测t_p时刻输出
    y(k) = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn);    %获取y(k)
    %if(k<71-First_t)
    %   y(k) = 0.7*y(k-1)+0.9*u(k-1)-0.4*u(k-2)+c*w(k);
    %else   %t=50时被控对象改变
    %   y(k) = -0.3*y(k-1)-0.7*y(k-2)+0.9*u(k-1)+0.7*u(k-2)+c*w(k);
    %end
    %------------------------------3.2 参数估计
    %----------1 更新Yd、Yabs_mtpb_2、Yr
    Yd = [y(k)-y(k-1);[Yd(1:na+N-1,1)]];
    Y = [y(k);[Y(1:na+N-1,1)]];
    Itmd_val = yr(k)*(1-soft_ele);
    Yr(1) = y(k)*soft_ele + Itmd_val;
    for j=2:N
        Yr(j) = Yr(j-1)*soft_ele + Itmd_val;
    end
    %----------2 计算φ
    phi = 0;
    for i=1:Nu
        phi = phi + Ud(N+1-i)*Q(i);
    end
    phi = phi * lambda;        %λ*Q*△u
    for i=N:-1:1
        phi = phi + P(i)*Y(N+1-i);
    end                        %P*y
    phi = phi - sum(P)*Y(N+1); %-sum(P)*y(k-N)
    %----------3 获取X(k-N),X(k)
    X_mns_N = [[Yd(N+1:N+na,1)];[Ud(N+1:N+nb,1)]];  %na+nb维向量
        %X(t-N)=[△y(k-N);...;△y(k-N-na+1);△u(t-1-N);...;△u(t-nb-N)]
    X = [[Yd(1:na,1)];[Ud(1:nb,1)]];
    %----------4 获取正则化因子
    regu_ele = max(abs(Y))*2;    %全称 regularization element 正则化因子
    Itmd_val = max(abs(U))*2;    %Itmd_val全称intermediate value 意为中间量
    regu_ele = max([regu_ele,Itmd_val,C]);
    %----------5 正则化
    phi = phi / regu_ele;
    X_mns_N = X_mns_N / regu_ele;
    %----------6 求ε(k),  f(M,ε(k))
    epsilon = phi - X_mns_N'*Theta - Ud(N)/regu_ele;
    if(epsilon>M)
        epsilon = epsilon - M;
    elseif(epsilon<-M)
        epsilon = epsilon + M;
    else
        epsilon = 0;
    end
    %----------7 更新θ(k)
    Theta = Theta + forget_ele*epsilon*X_mns_N/(1+X_mns_N'*X_mns_N);
    Theta_His(:,k) = Theta;
    %----------8 计算△u(k)
    ud = P*Yr - sum(P)*y(k) - X'*Theta;
    u(k) = u(k-1) + ud;  %求出u(t)
    if u(k)<u_min           %给u(t)一些限制条件
        u(k) = u_min;
        ud = u_min - u(k-1);
    elseif(u(k)>u_max)
        u(k) = u_max;
        ud = u_max - u(k-1);
    end
    ud_space(k)=ud;
    %----------9 更新Ud、U 
    Ud = [Ud(1:nb+N-1,1)];
    Ud = [ud;Ud];
    U = [U(1:nb+N-1,1)];
    U = [u(k);U];
end
run_time=toc;

figure(6);
uicontrol('Style','text','Position',[10 0 150 20],'String',['运行时间: ',num2str(run_time),'s']);
        %运行时间显示
subplot(211);
plot(t,yr,'b--',t,y,'r-');
xlim([First_t Last_t]);
%axis([First_t,Last_t,-10,60]);
legend('设定值yr','实际值y','Location','Best');
grid on;
title('基于被控对象前N个阶跃响应的GPC直接算法');
ylabel('输出y(k)');
subplot(212);
plot(t,ud_space);
xlim([First_t Last_t]);
grid on;
ylabel('控制增量u_d(k)');
xlabel('k');

%--------------------------------------------------5 编辑返回值
y = [y(1,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
u = [u(1,1-First_t:1+Last_t-First_t)];
ud_space = [ud_space(1,1-First_t:1+Last_t-First_t)];
Theta_His = [Theta_His(:,1-First_t:1+Last_t-First_t)];
DataTXT = fopen('.\GPC_simulate_data\GPC_d2__data.txt','w');%打开txt文件
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
