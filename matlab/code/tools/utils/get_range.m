function indices = get_range(idx, size)
% Helps convert cells of block matrices into giant matrix
% Returns the indices corresponding to the (idx)-th block,
% assuming blocks are uniformly of size (size)

indices = size*(idx-1)+1:size*(idx-1)+size;
