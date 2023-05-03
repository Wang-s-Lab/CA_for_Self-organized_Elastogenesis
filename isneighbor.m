%the following function determine whether a point in a matrix has at least
%a neighbor particle of group i.
function y=isneighbor(L,i)
    [nx,~]=size(L);
    n=nx;
    y=ones(n,n);
    j=1:n;
    lminus=[n,1:n-1];
    lplus=[2:n,1];
 y((L(lminus,j)~=i) + (L(lplus,j)~=i) + (L(j,lminus)~=i) + (L(j,lplus)~=i)==4)=0;
 y((L(lminus,j)~=i) + (L(lplus,j)~=i) + (L(j,lminus)~=i) + (L(j,lplus)~=i)~=4)=1;
end