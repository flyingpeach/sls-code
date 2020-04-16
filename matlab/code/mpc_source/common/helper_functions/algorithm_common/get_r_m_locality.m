function [r_loc, m_loc] = get_r_m_locality(sys, tFIR)

r_loc = cell(1, tFIR);
m_loc = cell(1, tFIR);

Comms_Adj = abs(sys.A)>0;
for t = 1:tFIR
    r_loc{t} = Comms_Adj^(d-1)>0;
    m_loc{t} = abs(sys.B2)'*r_loc{t}>0;
end

end