  %2019-03-08
  
%function：基本的GPC
%parameter: pN: the number of plants
%           PA: a group of parameters of the group of "y"(the output)
%               Plant1:a1,a2,a3,... Plant2:a1,a2,a3,... Plant3:a1,a2,a3,...
%           PB: a group of parameters of the group of "u"(the input)
%               Plant1:b0,b1,b2,... Plant2:b0,b1,b2,... Plant3:b0,b1,b2,...
%           Cn: the range of noise.（±Cn）
%return: y: 输出
%        u: 控制量
function [y,ud_space] = GPC_1(pN,PA,PB,na,nb,N1,Nu,lambda,soft_ele,Cn)
%############################## 用户编辑区 ##############################
[First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,N1);
%#######################################################################

%--------------------------------------------------1 表示其他已知量
u=zeros(1, Last_t+1-First_t);    %控制量
ud_space=zeros(1, Last_t+1-First_t);   %控制增量的存储
y=zeros(1,Last_t+1-First_t);    %输出量
w=rand(1,1+Last_t-First_t)*2-1;   %扰动。幅度范围（-1,1）的随机序列
for k = 2-First_t : 1+Last_t-First_t
    w(k) = w(k-1) + w(k);
end                             %表示ω(k)/?

A=[1,PA(1,:)];
B=PB(1,:);

%--------------------------------------------------2 离线计算控制率
[P,Alpha,Beta] = GPC_getCtrlRule(A,na,B,nb,N1,Nu,lambda);

%--------------------------------------------------3 控制过程
Ud = zeros(1,nb); Ud_p=1;  %存储△u(t-1),...,△u(t-nb)
Yr = zeros(N1,1);          %存储柔化后的设定值项
tic;
for k = 1-First_t : 1+Last_t-First_t   %时刻指针t, t_p=1-First_t对应实际时刻t=0
%for t = 0         :  Last_t
    %----------------------------------------3.1 检测k时刻输出
    y(k) = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn);    %获取y(k)
    %----------------------------------------3.2 更新控制量
    ud = 0;    %用于存储△u(t)
    %----- 设定值项
    mid = yr(k)*(1-soft_ele);
    Yr(1) = y(k)*soft_ele + mid;
    for j=2:N1
        Yr(j) = Yr(j-1)*soft_ele + mid;
    end
    ud = ud + P*Yr;    %其中P是横向量
    %for j=1:N1
        %if(t_p+j>Last_t+1-First_t)
        %    ud = ud + P(j)*yr(Last_t+1-First_t);
        %else
        %    ud = ud + P(j)*yr(t_p+j);
        %end
        %ud = ud + P(j)*yr(k);
    %end
    %----- 历史输出项
    for j=1:na+1
        ud = ud - Alpha(j)*y(k+1-j);
    end
    %----- 历史输入项
    for j=1:nb
        ud = ud - Beta(j)*Ud(mod(Ud_p-2+j,nb)+1);  %Ud_p指向△u(t-1),Ud_p+1指向△u(t-2)
    end    %已求出△u(t)
    %---------- 求出u(t)
    u(k) = u(k-1) + ud;
    if u(k)<u_min           %给u(t)一些限制条件
        u(k) = u_min;
        ud = u_min - u(k-1);
    elseif(u(k)>u_max)
        u(k) = u_max;
        ud = u_max - u(k-1);
    end
    ud_space(k)=ud;
    %---------- 更新△u(t-1)存储数组
    Ud_p = Ud_p-1;
    if Ud_p < 1
        Ud_p = nb;
    end
    Ud(Ud_p) = ud;     %以上5行。更新△u(t-1)
end
run_time=toc;

figure(3);
uicontrol('Style','text','Position',[10 0 150 20],'String',['运行时间: ',num2str(run_time),'s']);
        %运行时间显示
%set(gca,'position',[0.1 0.1 0.9 0.9]);
subplot(211);
plot(t,yr,'b--',t,y,'r-');
xlim([First_t Last_t]);
legend('设定值yr','实际值y','Location','Best');
grid on;
title('广义预测控制算法仿真（基本GPC，模型已知）');
ylabel('输出y(k)');
subplot(212);
plot(t,ud_space);
grid on;
xlim([First_t Last_t]);
ylabel('控制增量u_d(k)');
xlabel('k')

%--------------------------------------------------5 编辑返回值
y = [y(1,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
u = [u(1,1-First_t:1+Last_t-First_t)];
ud_space = [ud_space(1,1-First_t:1+Last_t-First_t)];
DataTXT = fopen('.\GPC_simulate_data\GPC_1__data.txt','w');%打开txt文件
fprintf(DataTXT,'%s\t%s\t%s\r\n','时刻k','输出y','控制增量ud');
data_len = length(y(:));
for j=1:data_len
    fprintf(DataTXT,'%d\t',j-1);
    fprintf(DataTXT,'%.3f\t%.3f\r\n',y(j),ud_space(j));
end
fclose(DataTXT);

%研究笔记：
%1、增大lambda，调节时间变长，但更平缓
%2、Nu越小，超调率越大，系统波动越大，响应时间有延时
