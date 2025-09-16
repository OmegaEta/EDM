function meanMatrix = meancellindouble(A)
    % Check if the input is empty.
    if isempty(A)
        error('The input cell array must not be empty.');
    end
    
    % Check the consistency of all matrix dimensions.
    try
        A_cat = cat(3, A{:});  % Convert the cell array into a three-dimensional array.
    catch
        error('The matrix dimensions are inconsistent and cannot be merged.');
    end
    
    %Calculate the mean along the third dimension
    meanMatrix = mean(A_cat, 3);
end