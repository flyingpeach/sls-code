function [r, s_r, c, s_c, r_loc, m_loc] = get_locality_constraints(sys, tFIR)

Nx = sys.Nx; Nu = sys.Nu;

r_loc = cell(1, tFIR);
m_loc = cell(1, tFIR);

Comms_Adj = abs(sys.A)>0;
for t = 1:tFIR
    r_loc{t} = Comms_Adj^(d-1)>0;
    m_loc{t} = abs(sys.B)'*r_loc{t}>0;
end

c     = cell(1, Nx);
s_c   = cell(1, Nx);

% Separate by columns (see columnwise_separability.m for details)
for i = 1:Nx
    c{i} = i;
    count = 0;
    for j = 1:tFIR+(tFIR-1)
        if j<=tFIR
            find_locR = find(r_loc{j}(:,i));
            for k =1:max(length(find_locR))
                count = count +1;
                s_c{i}(count) = find_locR(k)+(j-1)*Nx;
                if j == tFIR
                end
            end
        else
            find_locM = find(m_loc{j-tFIR}(:,i));
            for k =1:max(length(find_locM))
                count = count +1;
                s_c{i}(count) = find_locM(k)+(j-tFIR-1)*Nu+tFIR*Nx;
            end
        end
    end
end

r     = cell(1, Nx);
s_r   = cell(1, Nx);

% Separate by rows (see rowwise_separability.m for details)
k = 0;
for i = 1:Nx
    if mod(i, Nx/Nu) == 0 % Decide whether or not there is actuation
        s_r{i} = zeros(tFIR+(tFIR-1),Nx); % Prealocate the indices
        k = k+1;
        for j = 1:tFIR+(tFIR-1)
            if j<=tFIR
                r{i}(j) = Nx*(j-1) + i;
                s_r{i}(j,1:max(length(find(r_loc{j}(i,:))))) = find(r_loc{j}(i,:));
            else
                r{i}(j) = Nu*(j-tFIR-1) + Nx*tFIR + k;
                s_r{i}(j,1:max(length(find(m_loc{j-tFIR}(k,:))))) = find(m_loc{j-tFIR}(k,:));
            end
        end
    else
        s_r{i} = zeros(tFIR,Nx); % Prealocate the indices
        for j = 1:tFIR
            r{i}(j) = Nx*(j-1) + i;
            s_r{i}(j,1:max(length(find(r_loc{j}(i,:))))) = find(r_loc{j}(i,:));
        end
    end
    s_r{i}( :, ~any(s_r{i},1) ) = []; % Eliminate the columns with only zeros
end