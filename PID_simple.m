 %2019-04-02
 
 %function: ��ɢPID����
 %parameter: pN,PA,PB,na,nb: ���ض������
 %           Kp,Ti,Td: PID������Ki=Kp/Ti, Kd=Kp*Td
 function [y,ud_space]=PID_simple(pN,PA,PB,na,nb,Kp,Ti,Td,Cn)
%PA=[0.3 0.7],PB=[0.9 0.7]ʱ��Kp=0.2��Ti=1.65��Td=0.4���У��ȷ�Kp=0.43��Pc=3.3��
%############################## �û��༭�� ##############################
[First_t,Last_t,tim_len,t,yr,u_max,u_min] = GPC_ParameterEditSpace(pN,2);
%#######################################################################

%--------------------------------------------------1 ��ʾ������֪��
Ki = 0;
if(Ti<=65535)
    Ki = Kp/Ti;     %Ki
end
Kd = Kp*Td;         %Kd
u=zeros(1, Last_t+1-First_t);    %������
ud_space=zeros(1, Last_t+1-First_t);   %���������Ĵ洢
y=zeros(1,Last_t+1-First_t);    %�����
e=zeros(1,Last_t+1-First_t);    %ƫ������yr(k)-y(k)��
ei=0;                           %ƫ�����Ļ���ֵ
w=rand(1,1+Last_t-First_t)*2-1;   %�Ŷ������ȷ�Χ��-1,1�����������
for k = 2-First_t : 1+Last_t-First_t
    w(k) = w(k-1) + w(k);
end                             %��ʾ��(k)/?
%U=0;   %�洢��ǰʱ�̵Ŀ������� U=Kp*e(k)+Ki*��e(k)+Kd*��e(k) ��
%Ud=0;  %�洢��ǰʱ�̿������ı仯��( Ud=Kp*��e(k)+Ki*e(k)+Kd*��(��e(k)) )

%--------------------------------------------------3 ���ƹ���
tic;
for k = 1-First_t : 1+Last_t-First_t   %ʱ��ָ��t, t_p=1-First_t��Ӧʵ��ʱ��t=0
    %------------------------------3.1 ���t_pʱ�����
    y(k) = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn);    %��ȡy(k)
    %------------------------------3.2 ʵʱ���¿�����
    e(k) = yr(k) - y(k);    %e(k)��ʾƫ��
    %----- idea 1��λ����ʽ
    u(k) = Kp*e(k);                     %P
    ei = ei + e(k);u(k) = u(k) + ei*Ki; %I
    u(k) = u(k) + Kd*(e(k)-e(k-1));     %D
    %----- idea 2��������ʽ
    %Ud = Kp*(e(k)-e(k-1));
    %if(Ti<65536)  %I�����
    %    Ud = Ud + Kp/Ti*e(k);
    %end
    %if(Td>0)
    %    Ud = Ud + Kp*Td*(e(k)-2*e(k-1)+e(k-2));
    %end
    %u(k) = u(k-1) + Ud;
    if u(k)<u_min      %��u(t)һЩ��������
        u(k) = u_min;
    elseif(u(k)>u_max)
        u(k) = u_max;
    end
    ud_space(k)=u(k)-u(k-1);
end
run_time=toc;

figure(7);
uicontrol('Style','text','Position',[10 0 150 20],'String',['����ʱ��: ',num2str(run_time),'s']);
        %����ʱ����ʾ
subplot(211);
plot(t,yr,'b--',t,y,'r-');
legend('�趨ֵyr','ʵ��ֵy','Location','Best');
%axis([First_t Last_t -10 50]);
xlim([First_t Last_t]);   %ֻ���ƺ����귶Χ
grid on;
title('PID');
ylabel('���y(k)');
subplot(212);
plot(t,ud_space);
%axis([First_t Last_t -10 100]);
xlim([First_t Last_t]);
grid on;
ylabel('��������u_d(k)');
xlabel('k')

%--------------------------------------------------5 �༭����ֵ
y = [y(1,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
u = [u(1,1-First_t:1+Last_t-First_t)];
ud_space = [ud_space(1,1-First_t:1+Last_t-First_t)];
DataTXT = fopen('.\GPC_simulate_data\PID__data.txt','w');%��txt�ļ�
fprintf(DataTXT,'%s\t%s\t%s\r\n','ʱ��k','���y','��������ud');
data_len = length(y(:));
for j=1:data_len
    fprintf(DataTXT,'%d\t',j-1);
    fprintf(DataTXT,'%.3f\t%.3f\r\n',y(j),ud_space(j));
end
fclose(DataTXT);
