%Given matrix A outputs its matrix P of pivot values
function [P] = PivotValue(A)

[y,x]=size(A);%Records the size of A

P=zeros(y,x);%Outputs matrix of the correct size

%Computes value for columns
for k=1:x
    
    [Xc,Gc] = VecGCD(A(:,k).');%Computes gcd for current column

    if Gc%Checks the column was not a zero vector
    
		%If the first value of gcd scalar vector is zero changes it to an equivalent vector where the first entry is non-zero
        if Xc(1)==0
            temp=(A(1,k))/Gc;
            Xc=Xc*(temp+1);
            Xc(1)=-1;
        end;
		
		%Computes the values of the matrix if this column was the pivot
        ColVal=zeros(y,x);
        
        for i=1:y
            for j=1:x
                    ColVal(i,j)=abs(A(i,j)-((dot(Xc,A(:,j)))/(dot(Xc,A(:,k)))*A(i,k)));
            end;
        end;

        temp=max(max(ColVal));%Maximum value in column of the pivot matrix

        P(:,k)=P(:,k)+temp*ones(y,1);%Records the max value on the corresponding column of P
    
    else
    
        P(:,k)=P(:,k)+inf*ones(y,1);%Records zero column as infinite pivot value
        
    end;
    
end;

%Computes value for rows
for k=1:y
     
     [Xr,Gr] = VecGCD(A(k,:));

     if Gr%Checks the row was not a zero row
		
		%If the first value of gcd scalar vector is zero changes it to an equivalent vector where the first entry is non-zero
         if Xr(1)==0
             temp=(A(k,1))/Gr;
             Xr=Xr*(temp+1);
             Xr(1)=-1;
         end;
		
		%Computes the values of the matrix if this column were a pivot
         RowVal=zeros(y,x);
         for i=1:y
             for j=1:x
                 RowVal(i,j)=abs(A(i,j)-((dot(Xr,A(i,:)))/(dot(Xr,A(k,:)))*A(k,j)));
             end;
         end;

		%Checks the maximum value in the pivot matrix for this row
         temp=max(max(RowVal));

		%Multiplies the row of P by this max value
         for a=1:x
             P(k,a)=P(k,a)*temp;
         end;
     
     else
     
         for a=1:x
             P(k,a)=inf;%Records zero row as infinite pivot value 
         end;
         
     end;
         
 end;
