function value_on_node=project_cell2node_master(z,dz,value_in_cell)
    N=length(dz);
    value_on_node=zeros(N+1,1);
    value_on_node(1)=(value_in_cell(1)*(2*dz(1)+dz(2))-value_in_cell(2)*dz(1))/(dz(1)+dz(2));
    value_on_node(end)=(value_in_cell(end)*(2*dz(end)+dz(end-1))-value_in_cell(end-1)*dz(end))/(dz(end)+dz(end-1));
    
    for i=2:N
        value_on_node(i)=(value_in_cell(i-1)*dz(i)+value_in_cell(i)*dz(i-1))/(dz(i-1)+dz(i));
    end
end