% Function to make sure that folder actually exists.  Warn user if it doesn't.

function f = foldCheck(fname)
if ~isdir(fname)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', fname);
    uiwait(warndlg(errorMessage));
    return;
end

end
