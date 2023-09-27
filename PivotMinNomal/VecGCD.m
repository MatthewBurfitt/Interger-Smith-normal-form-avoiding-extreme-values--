%Given a vector V outputs gcd G and vector of scales X whose scalar product with V is G
function [X,G] = VecGCD(V)

s=size(V,2);%number of elements in V

P=eye(s);
%For recording intermediary values of X
minV=[1,inf];

temp=0;

neg=zeros(1,s);%for recording sign changes

%ensures V is non-negative integer vector and vectors where the sign changes
for i=1:s
    if V(i)<0
        V(i)=-V(i);
        neg(i)=1;
    end;
end;

%Checks for exceptional vase where V is the zero vector   
if V==zeros(1,s)
    X=zeros(1,s);
    G=0;
else
%Computes G and X with the Euclidean algorithm
while minV(1)

    minV=[0,inf];

	%Finds the smallest value i in V
    for i=1:s
        if V(i)
            if V(i)<minV(2)
                minV=[i,V(i)];
            end;
        end;
    end;
    
	%If the minimum positive value is unchanged this is the gcd and the procedure terminates
    if temp==minV(2)
        G=minV(2);
        X=P(minV(1),:);
        break
    end;
    
	%reduce the vector V modulo its minimum value and record what was done in P
    if minV(1)
        for i=1:s
            if i==minV(1)
            else
                f=floor(V(i)/minV(2));
                V(i)=V(i)-f*minV(2);
                P(i,:)=P(i,:)-f*P(minV(1),:);
            end;
        end;
    end;
    
    temp=minV(2);

end;
end;

%assigns the correct sign to elements of X
for i=1:s
    X(i)=X(i)*(-1)^(neg(i));
end;


