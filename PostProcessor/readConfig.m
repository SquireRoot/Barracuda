function [out] = readConfig()

fh = fopen('../DATA/config.txt');
file = textscan(fh, '%s %s %d', 6);
fclose(fh);
data = file{3};

out.sizeX = data(1);
out.sizeY = data(2);
out.sizeZ = data(3);
out.iterations = data(4);
out.outputInterval = data(5);
out.outputBegin = data(6);

end

