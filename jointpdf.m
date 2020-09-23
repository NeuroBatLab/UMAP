function [X,Y,pdf_normalize]=jointpdf(data_x,data_y,N,plon)

stepx=(max(data_x)-min(data_x))/N;
stepy=(max(data_y)-min(data_y))/N;

x = min(data_x):stepx:max(data_x);
y = min(data_y):stepy:max(data_y); 

[X,Y] = meshgrid(x,y); 

pdf = hist3([data_x , data_y],{x y});
pdf_normalize = (pdf'./length(data_x)); 
   if plon                                     
figure()
surf(X,Y,pdf_normalize)
   end
end
