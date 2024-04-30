%2019-04-16

%function：输入被控对象参数以及已知量，求得y(k)
%parameter: 参数含义与GPC各算法中的含义均相同
function yk = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn)
plt = fix((k+First_t-1)/tim_len)+1;  %plt表示被控对象序号
if(plt>3)
    plt=3;   %当需要仿真多个被控对象时，第4、5...个被控对象参数只能与第3个相同
end
len_a=length(PA(1,:));
len_b=length(PB(1,:));
yk = Cn*w(k);
for i=1:len_a
    yk = yk - PA(plt,i)*y(k-i);
end
for i=1:len_b
    yk = yk + PB(plt,i)*u(k-i);
end    %获取y(k)
%if(k+First_t-1==240)   %手动添加脉冲干扰
%    yk=yk+10;
%end