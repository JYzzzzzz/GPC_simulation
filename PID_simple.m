 %2019-04-02
 
 %function: 离散PID仿真
 %parameter: pN,PA,PB,na,nb: 被控对象参数
 %           Kp,Ti,Td: PID参数，Ki=Kp/Ti, Kd=Kp*Td
 function [y,ud_space]=PID_simple(pN,PA,PB,na,nb,Kp,Ti,Td,Cn)
%PA=[0.3 0.7],PB=[0.9 0.7]时，Kp=0.2，Ti=1.65，Td=0.4可行（等幅Kp=0.43，Pc=3.3）
%############################## 用户编辑区 ##############################
[First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,2);
%#######################################################################

%--------------------------------------------------1 表示其他已知量
Ki = 0;
if(Ti<=65535)
    Ki = Kp/Ti;     %Ki
end
Kd = Kp*Td;         %Kd
u=zeros(1, Last_t+1-First_t);    %控制量
ud_space=zeros(1, Last_t+1-First_t);   %控制增量的存储
y=zeros(1,Last_t+1-First_t);    %输出量
e=zeros(1,Last_t+1-First_t);    %偏差量（yr(k)-y(k)）
ei=0;                           %偏差量的积分值
w=rand(1,1+Last_t-First_t)*2-1;   %扰动。幅度范围（-1,1）的随机序列
for k = 2-First_t : 1+Last_t-First_t
    w(k) = w(k-1) + w(k);
end                             %表示ω(k)/?
%U=0;   %存储当前时刻的控制量（ U=Kp*e(k)+Ki*∑e(k)+Kd*△e(k) ）
%Ud=0;  %存储当前时刻控制量的变化量( Ud=Kp*△e(k)+Ki*e(k)+Kd*△(△e(k)) )

%--------------------------------------------------3 控制过程
tic;
for k = 1-First_t : 1+Last_t-First_t   %时刻指针t, t_p=1-First_t对应实际时刻t=0
    %------------------------------3.1 检测t_p时刻输出
    y(k) = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn);    %获取y(k)
    %------------------------------3.2 实时更新控制率
    e(k) = yr(k) - y(k);    %e(k)表示偏差
    %----- idea 1：位置形式
    u(k) = Kp*e(k);                     %P
    ei = ei + e(k);u(k) = u(k) + ei*Ki; %I
    u(k) = u(k) + Kd*(e(k)-e(k-1));     %D
    %----- idea 2：增量形式
    %Ud = Kp*(e(k)-e(k-1));
    %if(Ti<65536)  %I项存在
    %    Ud = Ud + Kp/Ti*e(k);
    %end
    %if(Td>0)
    %    Ud = Ud + Kp*Td*(e(k)-2*e(k-1)+e(k-2));
    %end
    %u(k) = u(k-1) + Ud;
    if u(k)<u_min      %给u(t)一些限制条件
        u(k) = u_min;
    elseif(u(k)>u_max)
        u(k) = u_max;
    end
    ud_space(k)=u(k)-u(k-1);
end
run_time=toc;

figure(7);
uicontrol('Style','text','Position',[10 0 150 20],'String',['运行时间: ',num2str(run_time),'s']);
        %运行时间显示
subplot(211);
plot(t,yr,'b--',t,y,'r-');
legend('设定值yr','实际值y','Location','Best');
%axis([First_t Last_t -10 50]);
xlim([First_t Last_t]);   %只限制横坐标范围
grid on;
title('PID');
ylabel('输出y(k)');
subplot(212);
plot(t,ud_space);
%axis([First_t Last_t -10 100]);
xlim([First_t Last_t]);
grid on;
ylabel('控制增量u_d(k)');
xlabel('k')

%--------------------------------------------------5 编辑返回值
y = [y(1,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
u = [u(1,1-First_t:1+Last_t-First_t)];
ud_space = [ud_space(1,1-First_t:1+Last_t-First_t)];
DataTXT = fopen('.\GPC_simulate_data\PID__data.txt','w');%打开txt文件
fprintf(DataTXT,'%s\t%s\t%s\r\n','时刻k','输出y','控制增量ud');
data_len = length(y(:));
for j=1:data_len
    fprintf(DataTXT,'%d\t',j-1);
    fprintf(DataTXT,'%.3f\t%.3f\r\n',y(j),ud_space(j));
end
fclose(DataTXT);
