%Takes matrix A and put it in its Smith normal form over filed of prime characteristic p, outputting the number of entries on the leading diagonal
function [U] = ModuloNormalForm(A,p)

U=0;

h=size(A,1);%Height of A
w=size(A,2);%Width of A

L=min(h,w);%The size of the leading diagonal

%The normal form procedure moves the current position along the leading diagonal
for a=1:L
    
	%Reduces the matrix to its simplest integer representative modulo p
    for i=a:h
        for j=a:w
            if A(i,j)>0
                A(i,j)=A(i,j)-floor(A(i,j)/p)*p;
            else
                A(i,j)=A(i,j)-floor(A(i,j)/p)*p;
            end;
        end;
    end;
    
    done=1;
    
	%Checks to see if the current row and column are zero an if so proceeds to the next position on the leading diagonal
    if A(a,:)==zeros(1,w)
        if A(:,a)==zeros(h,1)
            done=0;
        end;
    end;
    
%Uses integral row and column operations to reduce the current position to the greatest common divisor of its row, then all other entries to zero
    while done
        
		%Moves the smallest positive integer in the current row or column to the current position        
        Low=[A(a,a),a,0];
                
        if Low(1)
        else
             Low(1)=inf;
        end;
                
        for i=a+1:h
            if A(i,a)
                if A(i,a)<Low(1)
                    Low=[A(i,a),i,0];
                end;
            end;
        end;
                
        for i=a+1:w
            if A(a,i)
                if A(a,i)<Low(1)
                    Low=[A(a,i),i,1];
                end;
            end;
        end;
                
        if Low(3)
            A(:,[a,Low(2)])=A(:,[Low(2),a]);
        else
            A([a,Low(2)],:)=A([Low(2),a],:);
        end;
                
        done=0;
         
		%Reduces all non-zero entries in the current column by the integer in the current position
        for i=a+1:h
            if A(i,a)
                A(i,:)=A(i,:)-floor(A(i,a)/A(a,a))*A(a,:);
            end;
            if A(i,a)
                done=1;
            end;
        end;
        
		%Reduces all non-zero entries in the current row by the integer in the current position		
        for i=a+1:w
            if A(a,i)
                A(:,i)=A(:,i)-floor(A(a,i)/A(a,a))*A(:,a);
            end;
            if A(a,i)
                done=1;
            end;
        end;
                
    end;
	
	%If no reduction took place than move the current position on the leading diagonal otherwise repeat from finding the smallest entry
        
end;

%Reduce the final diagonal from the matrix modulo p
for i=a:h
    for j=a:w
        if A(i,j)>0
            A(i,j)=A(i,j)-floor(A(i,j)/p)*p;
        else
            A(i,j)=A(i,j)+floor(A(i,j)/p)*p;
        end;
    end;
end;

temp=0;

%Counts the number of non-zero entries on the leading diagonal of the normal form matrix
for i=1:L
    if A(i,i)
        temp=temp+1;
    end;
end;

U=w-temp;%Outputs the number of non-zero entries on the leading diagonal of the normal form matrix
