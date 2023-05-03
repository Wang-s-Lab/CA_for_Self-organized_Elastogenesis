function [output]=randommove(cells,movestatusL,coNei)
    [n1,n2]=size(cells);
    n=n1;
    BW=cells;
    [B,L] = bwboundaries(BW,coNei,'noholes');
    sums=zeros(n1,n2);
    for rightbno=1:n
        if L(rightbno,n)~=0 && L(rightbno,1)~=0 && L(rightbno,n)~=L(rightbno,1)
            Ltemp=L(rightbno,n);
            L(L==Ltemp)=L(rightbno,1);
        end
    end
    
    for bottombno=1:n
        if L(n,bottombno)~=0 && L(1,bottombno)~=0 && L(n,bottombno)~=L(1,bottombno)
            Ltemp=L(n,bottombno);
            L(L==Ltemp)=L(1,bottombno);
        end
    end
    
    for i=1:max(L(:))
        if movestatusL(i)==1
            randsteppara=randi([1,4]);
            location=find(L==i);
            [R,C]=ind2sub(size(L),location);
            for j=1:length(R)
                cells(R(j),C(j))=0;
            end
            
            for j=1:length(R)
                if randsteppara==1
                    xcoordinate=R(j)+1;
                    ycoordinate=C(j);
                end
                if randsteppara==2
                    xcoordinate=R(j)-1;
                    ycoordinate=C(j);
                end
                if randsteppara==3
                    xcoordinate=R(j);
                    ycoordinate=C(j)+1;
                end
                if randsteppara==4
                    xcoordinate=R(j);
                    ycoordinate=C(j)-1;
                end
                
                %periodic BCs
                for coordinateindex=1:length(xcoordinate)
                    if xcoordinate(coordinateindex)>n
                        xcoordinate(coordinateindex)=xcoordinate(coordinateindex)-n;
                    end
                    if xcoordinate(coordinateindex)<1
                        xcoordinate(coordinateindex)=xcoordinate(coordinateindex)+n;
                    end
                    if ycoordinate(coordinateindex)>n
                        ycoordinate(coordinateindex)=ycoordinate(coordinateindex)-n;
                    end
                    if ycoordinate(coordinateindex)<1
                        ycoordinate(coordinateindex)=ycoordinate(coordinateindex)+n;
                    end
                end
                cells(xcoordinate,ycoordinate)=1;
            end
        elseif movestatusL(i)==0
        end
    end
    output=cells;
end