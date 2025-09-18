function Poles = getPoles(hm)
Poles = zeros(2*hm,1);
for index = 1:hm
    theta = 2*(pi*index)^0.5;
    r = (index/(hm+1))^0.5;
    SpiralPole1 = r*exp(1i*theta);
    SpiralPole2 = r*exp(-1i*theta);
    Poles(2*index-1) = (SpiralPole1-1)/(SpiralPole1+1);
    Poles(2*index) = (SpiralPole2-1)/(SpiralPole2+1);
end
end