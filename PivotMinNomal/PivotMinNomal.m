%Given a matrix A finds its Smith Normal form in a way that attempts to minimize the magnitude of intermediary values
function [A] = PivotMinNomal(A)

[y,x]=size(A);%Records the size of A

max=min(x,y);%size of the leading diagonal

for i=1:max
%i
    null=1;
    
    for a=i:x %check to see if all remaining entries are zero
        for b=i:y
            if A(b,a)
                null=0;
                break
            end;
        end;
        if null
        else
            break
        end;
    end;
    
    if null
        break
    end;
    
	%takes B the part of the matrix which we still need to reduce
	
    B=zeros(y-i+1,x-i+1);
    
    for a=i:x
        for b=i:y
            B(b-i+1,a-i+1)=A(b,a);
        end;
    end;
       
    B=PivotValue(B);
    
    MinPiv=[1,1,inf];
    
    for a=1:size(B,1) %finds none-zero value with smallest pivot value
        for b=1:size(B,2)
            if A(i+a-1,i+b-1)
                if B(a,b)<MinPiv(3)
                    MinPiv=[a,b,B(a,b)];
                end;
            end;
        end;
    end;
    
    p=MinPiv(1);
    q=MinPiv(2);
    
    A(:,[i,q+i-1])=A(:,[q+i-1,i]);
    A([i,p+i-1],:)=A([p+i-1,i],:);
    
    %now perform GCD reduction on the first row column for the top left position.
    in=1;
    
    while in
       
        
        if A(i,i)<0
            A(i,:)=-1*A(i,:);
        end;

        for a=i+1:y
            if A(a,i)<0
                A(a,:)=-1*A(a,:);
            end;
            A(a,:)=A(a,:)-floor(A(a,i)/A(i,i))*A(i,:);
        end;

        for a=i+1:x
            if A(i,a)<0
                A(:,a)=-1*A(:,a);
            end;
            A(:,a)=A(:,a)-floor(A(i,a)/A(i,i))*A(:,i);
        end;
        
        %check to see if all first row and column are zero except top left.
        out=1;
        
        for a=i+1:y
            if A(a,i)
                out=0;
            end;
        end;
        
        for a=i+1:x
            if A(i,a)
                out=0;
            end;
        end;
        
        if out
           break 
        end;
        
        %finds new pivot in fist row or column and repeat reduction
        
        B=zeros(y-i+1,x-i+1);
    
        for a=i:x
            for b=i:y
                B(b-i+1,a-i+1)=A(b,a);
            end;
        end;
        
        B=PivotValue(B);
        
        V=B(:,1).';
        
        H=B(1,:);
        
        piv=[1,1,inf];
        
        p=abs(A(i,:));
        q=abs(A(:,i)).';
        
        P=p(1);
        
        for a=2:size(p,2)
            if P<p(a)
                P=p(a);
            end;
        end;
        
        for a=1:size(q,2)
            if P<q(a)
                P=q(a);
            end;
        end
        
        temp=P;
        
        U=0;
        
        for a=i:x
            if abs(A(i,a))==temp
                U=U+1;
            end;
        end;
        
        for a=i+1:y
            if abs(A(a,i))==temp
                U=U+1;
            end;
        end;
        
        for a=1:size(V,2) %find lowest pivot value
            if A(a+i-1,i)
                if V(a)<piv(3)
                    if temp>abs(A(a+i-1,i))
                        piv=[1,a,V(a)];
                    else
                        if U>1
                            piv=[1,a,V(a)];
                        end;
                    end;
                end;
            end;
        end;
    
        for a=1:size(H,2)
            if A(i,a+i-1)
                if H(a)<piv(3)
                    if temp>abs(A(i,a+i-1))
                        piv=[0,a,H(a)];
                    else
                        if U>1
                            piv=[0,a,H(a)];
                        end;
                    end;
                end;
            end;
        end;
   
        if piv(1)
            A([i,piv(2)+i-1],:)=A([piv(2)+i-1,i],:);
        else
            A(:,[i,piv(2)+i-1])=A(:,[piv(2)+i-1,i]);
        end;
        
    end;
    
end;

%rearranges elements on diagonal smallest towards top left.
swap=1;

while swap
    swap=0;
    for i=1:max-1
        if A(i+1,i+1)==0
            break
        end;
        if A(i,i)>A(i+1,i+1)
            temp=A(i,i);
            A(i,i)=A(i+1,i+1);
            A(i+1,i+1)=temp;
            swap=1;
        end;
    end;
end;

