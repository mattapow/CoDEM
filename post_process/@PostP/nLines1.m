function N = nLines1(fp)
% returns the number of newline characters (ASCII 10) in the file

fID = fopen(fp);
try
    s = fread(fID);
catch
    error('Error reading filename: %s', fp)
end
N = sum(s==10);
fclose(fID);

end
