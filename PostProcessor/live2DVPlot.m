model = readConfig();

makeVideo = false;
zLayer = 1;

if makeVideo
    vh = VideoWriter('TRTP.avi');
    vh.FrameRate = 5;
    open(vh);
end   

%fig = figure;
for k = model.outputBegin:model.outputInterval:model.iterations
  filename = sprintf('../DATA/v/file%07d.dat', k);
  fh = fopen(filename);
  while fh == -1
      pause(1);
      fh = fopen(filename);
  end
  
  data = fread(fh, model.sizeX.*model.sizeY.*model.sizeZ.*3, 'float32');
  x = 1:model.sizeX;
  y = 1:model.sizeY;
  z = 1:model.sizeZ;
  data = reshape(data, [3, model.sizeX, model.sizeY, model.sizeZ]);
  xData = flipud(squeeze(data(1, :, :, zLayer))');
  yData = flipud(squeeze(data(2, :, :, zLayer))');
  zData = flipud(squeeze(data(3, :, :, zLayer))');
  
  v = sqrt(xData.^2 + yData.^2 + zData.^2);
  
  clf;
  imagesc(v, [0, 0.04]);
  colorbar;
  %hold on
  %y = 200:2:300;
  %x = ones(1, length(y));
  %streamline(xData,yData,x,y);
  %quiver(xData, yData);
  titlename = sprintf('Velocity at time step %d',k) ;
  title(titlename)
  axis equal
  set(gcf, 'Position',  [100, 100, 612, 512])
  
  pause(0.01);
  
  if makeVideo
      frame = getframe(fig);
      writeVideo(vh, frame);
  end
  
end

if makeVideo
    close(vh);
end