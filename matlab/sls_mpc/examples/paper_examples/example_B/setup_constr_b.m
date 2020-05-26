function K = setup_constr_b(Nx)

for j = 1:2:2*(Nx-1)
    K(j,j)     = 1; 
    K(j,j+2)   = -1;
    K(j+1,j)   = -1; 
    K(j+1,j+2) = 1;
end
K            = K(1:Nx,1:Nx); 
K(Nx-1:Nx,:) = zeros(2,Nx);

end