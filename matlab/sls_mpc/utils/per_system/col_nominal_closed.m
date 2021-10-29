function psi = col_nominal_closed(phi, lamb, zab, myeye, zabi)
% phi   : column of phi
% lamb  : column of Lambda
% zab   : column of ZAB (part of sls constraints)
% myeye : column of shifted identity matrix (part of sls constraints)
% zabi  : pseudoinverse of zab

psi = (phi + lamb) + zabi*(myeye - zab*(phi + lamb));

end