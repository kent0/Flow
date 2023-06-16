function res=lag(alag,a);
    res=alag;
    for i=0:size(res,2)-2
        res(:,end-i)=res(:,end-i-1);
    end
    res(:,1)=a;
end
