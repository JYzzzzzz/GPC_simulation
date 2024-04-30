%2019-04-16

%function�����뱻�ض�������Լ���֪�������y(k)
%parameter: ����������GPC���㷨�еĺ������ͬ
function yk = GPC_get_yk(First_t,k,tim_len,u,y,w,PA,PB,Cn)
plt = fix((k+First_t-1)/tim_len)+1;  %plt��ʾ���ض������
if(plt>3)
    plt=3;   %����Ҫ���������ض���ʱ����4��5...�����ض������ֻ�����3����ͬ
end
len_a=length(PA(1,:));
len_b=length(PB(1,:));
yk = Cn*w(k);
for i=1:len_a
    yk = yk - PA(plt,i)*y(k-i);
end
for i=1:len_b
    yk = yk + PB(plt,i)*u(k-i);
end    %��ȡy(k)
%if(k+First_t-1==240)   %�ֶ�����������
%    yk=yk+10;
%end