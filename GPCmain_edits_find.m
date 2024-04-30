%2019-04-17

%function：进入GPC仿真主界面时，辅助编辑框中数据的初始化、替换操作
%parameter: retType: 返回值的形式。
%                    dat: EditData中的数据，或
%                    idx: EditData对应的下标
function retData = GPCmain_edits_find(EditName,EditData,str,retType)
EditLen = length(EditName(:));
retData = 0;
for i=1:EditLen  %顺序遍历
    if(strcmp(EditName(i),str)==1) %找到对应的一项
        if(strcmp(retType,'dat')==1) %返回EditData中的数据
            retData = EditData(i);
        elseif(strcmp(retType,'idx')==1) %返回EditData对应的下标
            retData = i;
        end
        return; %提前跳出程序
    end
end