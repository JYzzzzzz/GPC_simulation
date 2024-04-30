  %2018-03-08
%function：绘制被控对象的阶跃响应
%parameter: 
function GPC_getStepResponse(pN,PA,PB)
%---------- 用户编辑区 ----------
First_t = -10;  %开始时刻
Last_t = 40;    %结束时刻
if(pN>3)        %pN上限
    pN=3;
end
%===============================
t = First_t : Last_t;
u=stepfun(t,0);           %单位阶跃函数
y=zeros(1,Last_t+1-First_t); %阶跃响应存储位置
y_out=zeros(pN,Last_t+1-First_t);

figure(2);
for i=1:pN
    for k = 1-First_t:Last_t+1-First_t   %时刻指针, k=1-First_t对应t=0
        y(k) = -PA(i,1)*y(k-1)-PA(i,2)*y(k-2)-PA(i,3)*y(k-3)+PB(i,1)*u(k-1)+PB(i,2)*u(k-2)+PB(i,3)*u(k-3);
        y_out(i,k) = y(k);
    end
    subplot(pN,1,i);
    plot(t,y);
    grid on;
    tmpStr='Plant 1';
    tmpStr(7)=num2str(i);
    title(tmpStr)
end

y_out = [y_out(1:pN,1-First_t:1+Last_t-First_t)];   %k = 0 ~ Last_t
DataTXT = fopen('.\GPC_simulate_data\StepResponse__data.txt','w');%打开txt文件
fprintf(DataTXT,'%s\t%s\r\n','时刻k','被控对象阶跃响应');
data_len = Last_t+1;
for j=1:data_len
    fprintf(DataTXT,'%d\t',j-1);
    for i = 1:pN
        fprintf(DataTXT,'%.4f\t',y_out(i,j));
    end
    fprintf(DataTXT,'\r\n');
end
fclose(DataTXT);
