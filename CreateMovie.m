clear all;
close all;
clc;
format compact;

writerObj = VideoWriter('Movie.mp4');
folder='testje'
List = dir(folder);
open(writerObj);
for K = 1:1:size(List,1)
  K
  Name = List(K).name;
  if (strfind(Name, '.png'))
    filename = strcat(folder, '/', Name);
    thisimage = imread(filename);
    writeVideo(writerObj, thisimage);
  end
end
close(writerObj);