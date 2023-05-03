function [a,b]=boundary(a,b,n)
for i=1:length(a)
   if a(i)>n
       a(i)=a(i)-n;
   elseif a(i)<1
       a(i)=a(i)+n;
   end
   
   if b(i)>n
       b(i)=b(i)-n;
   elseif b(i)<1
       b(i)=b(i)+n;
   end
       
end
 
end