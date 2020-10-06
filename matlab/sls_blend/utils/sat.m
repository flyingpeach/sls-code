function val_sat = sat(val, eta)
% clamps val to [-alpha, eta]

val_sat = val;

for i=1:length(val)
    if val(i) > eta
        val_sat(i) = eta;
    elseif val(i) < -eta
        val_sat(i) = -eta; 
    end
end

end