P=p';
N=length(P);  % N = total no. of Particles
B=25;L=70;  % Only for rectangular Lattice when bottom left corner is at origin, and L & B are length and breadth respectively
Position(N,2)=0;
maxNeighbour=8;
Neighbour(N,8)=0;
Neighbours_more_than_6=0;
Neighbours_less_than_6=0;
EqDistance(N,8)=0;
T=t';
T=[T(:,1) T(:,2) T(:,3)];
size_T=size(T);
for i=1:1:N
    y=mod(find(T==i),size_T(1));
    for j=1:1:length(y)
        if y(j)==0
            y(j)=size_T(1);
        end
    end
    temp_Neighbour=union(T(y,:),T(y,:));
    temp_Neighbour(temp_Neighbour==i)=[];
    if length(temp_Neighbour)>8
        fprintf('%i %s \n',i,'Neighbour of 1 particle exceeding limit of 8')
        maxNeighbour=max(length(temp_Neighbour),maxNeighbour);
    end
    if length(temp_Neighbour)>6
        %fprintf('%i %s \n',i,'Neighbour of 1 particle exceeding limit of 8')
        Neighbours_more_than_6=Neighbours_more_than_6+1;
    end
    if length(temp_Neighbour)<6
        %fprintf('%i %s \n',i,'Neighbour of 1 particle exceeding limit of 8')
        Neighbours_less_than_6=Neighbours_less_than_6+1;
    end
    
    Neighbour(i,1:length(temp_Neighbour))=temp_Neighbour;
        
    Position(i,:)=P(i,:);
end

% write position of all the particles to a file
% write Neighbours of all the particles to a file
file_Pos=fopen('position_matlab.dat','w');
file_Neighbour=fopen('neighbour_matlab.dat','w');
for i=1:1:N
fprintf(file_Pos,'%.10e \t %.10e \r\n',Position(i,:));
    for j=1:1:maxNeighbour
    fprintf(file_Neighbour,'%i \t',Neighbour(i,j));
    end
    fprintf(file_Neighbour,'\r\n');
end
fclose(file_Pos);
fclose(file_Neighbour);

% To calculate equillibrium distances
 file_EqDist=fopen('EqDistance_matlab.dat','w');
for i=1:1:N
    for j=1:1:maxNeighbour
         if Neighbour(i,j)~=0
            x= P(Neighbour(i,j),1)-P(i,1); 
            y= P(Neighbour(i,j),2)-P(i,2);
            EqDistance(i,j)=sqrt(x*x+y*y);
            fprintf(file_EqDist,'%.10e \t',EqDistance(i,j));
         end
    end
    fprintf(file_EqDist,'\r\n');
 end
fclose(file_EqDist);

%To calculate Particles at top and lower boundary
file_Top=fopen('topBoundary_matlab.dat','w');
file_Bottom=fopen('bottomBoundary_matlab.dat','w');
topBoundary=find(P(:,2)==B);
bottomBoundary=find(P(:,2)==0);
for i=1:1:length(topBoundary)
fprintf(file_Top,'%i \n',topBoundary(i));
end
for i=1:1:length(bottomBoundary)
fprintf(file_Bottom,'%i \n',bottomBoundary(i));
 end
fclose(file_Top);
fclose(file_Bottom);
