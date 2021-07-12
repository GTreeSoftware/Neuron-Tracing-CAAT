function [ Paths,new_branch_points,new_end_points,new_con_num,new_con_nodes ] = Build_Tracing_Graph( connect_links,connect_num,connect_nodes,line_threshold )
% this graph is used for tracing, start from one point and not to the end
% branch_points=find(connect_num>2);
end_points=find(connect_num==1);

% record the visit points
visit_connect_num=zeros(1,length(connect_num));
Paths=cell(length(end_points),1);
tmp_ep=end_points;

% record whether the start point generate multiple paths
ini_pnum=zeros(1,length(connect_num));

% begin to find the path for all
path_num=0;
while ~isempty(tmp_ep)
    % record path lenth
    p_length=0;
    
    % for the first time: from the last points 
    path_record=zeros(length(connect_num),1);
    start_points=tmp_ep(end);
    
    % mark the result
    visit_connect_num(start_points)=visit_connect_num(start_points)+1;
    tmp_ep(end)=[];
    
    % record the first path
    path_record(p_length+1)=start_points;
    p_length=p_length+1;
    
    cur_point=start_points;
    % begin to search
    next_points=connect_links(start_points);
    
    while 1
        % judge whether the points has been visited and not branch points
        if next_points==0
            break;      
        end
        
        if next_points==1
            % the start point
            ini_pnum(path_num+1)=2;    
        end
        
        if (visit_connect_num(next_points)==1)&&(connect_num(next_points)<3)
            break;           
        end
        
        % last point is branch point and current point is branch point
        % after the first line is searched 
        if path_num>=1
            if (visit_connect_num(next_points)>=1)&&(visit_connect_num(cur_point)>1)&&(connect_num(next_points)>=3)&&(connect_num(cur_point)>=3)
                break;           
            end
        end
        
        % judge whether the point is end points or not
        eflag=connect_num(next_points);
        
        % save this result
        path_record(p_length+1)=next_points;
        p_length=p_length+1;
        visit_connect_num(next_points)=visit_connect_num(next_points)+1;
        
        cur_point=next_points;
        
        % end points
        if eflag==1           
            % remove the points from the endpoints
            ind1=find(tmp_ep==next_points);
            tmp_ep(ind1)=[];
            
            % jump out of this iteration
            break;
        else
            % begin to search for the new points
             next_points=connect_links(next_points);    
        end
         
    end
    
    % judge the length of path and decide whether to remain or not
    path_record=path_record(1:p_length);
    
    % prunning short branches for simplify
    if p_length>line_threshold
        path_num=path_num+1;
        % fliplr the path
        path_record=flipud(path_record);
        Paths{path_num,1}=path_record;
        
        % decide the line from the first is kept
        if ini_pnum(path_num)>0
            ini_pnum(path_num)=1;
        end
        
    end
end
Paths=Paths(1:path_num,1);

% judge whether the seed point generate multiple branches
ini_ind1=find(ini_pnum==1);
ini_len1=length(ini_ind1);
if ini_len1>1
    tmp1=Paths{ini_ind1(1)};
    tmp2=Paths{ini_ind1(end)};
    tmp_all=[flipud(tmp2(2:end));tmp1];
    Paths{ini_ind1(1),1}=tmp_all;
    Paths(ini_ind1(2))=[];
    path_num=path_num-1;
end

% renew the branch points and end points
points_num=length(connect_num);
new_con_num=zeros(points_num,1);
new_con_nodes=zeros(size(connect_nodes));

% have branch points longer than the defined points because of two branch
% points are close and correspongding to the same lines
for i=1:path_num
    path_tmp=Paths{i,1};
    for j=2:size(path_tmp,1)
        lind1=path_tmp(j-1);
        cind1=path_tmp(j);
        
        % record the connection
        new_con_num(lind1)=new_con_num(lind1)+1;
        new_con_num(cind1)=new_con_num(cind1)+1;
        new_con_nodes(lind1,new_con_num(lind1))=cind1;
        new_con_nodes(cind1,new_con_num(cind1))=lind1;
                
    end
end

new_con_nodes=new_con_nodes(:,1:max(new_con_num));
new_branch_points=find(new_con_num>2);
new_end_points=find(new_con_num==1);


end

