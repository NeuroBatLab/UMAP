function plot_spec_cluster(DATA,bin,fhigh,compN,clN,mi,rep)

tits={'training';'test'};
%revert each line back into spectrogram
if mi==1
figure('Name',['Spec: Bin' num2str(bin) '-' tits{mi} '-Clus#' num2str(clN) '-Rep#' num2str(rep)],'NumberTitle','off');
else
 figure('Name',['Spec: Bin' num2str(bin) '-' tits{mi} '-Comp#' num2str(compN) '-Clus#' num2str(clN) '-Rep#' num2str(rep)],'NumberTitle','off');   
end
reduce=1;
%plot only subset if very large cluster
if size(DATA,1)>150 && size(DATA,1)<=300
    reduce=2;
elseif size(DATA,1)>300 && size(DATA,1)<=450
    reduce=4;
elseif size(DATA,1)>450 && size(DATA,1)<=600
    reduce=10;
elseif size(DATA,1)>600
    reduce=50;
end
for i=1:2:floor(size(DATA,1)/reduce)
    vec=reshape(DATA(i,:),length(DATA(i,:))/bin,bin);
    to=linspace(0,bin,size(vec,2));
    fo=linspace(0,fhigh,size(vec,1));
    subplot(ceil((size(DATA,1)/reduce)/(ceil((size(DATA,1)/reduce)/2))),ceil((size(DATA,1)/reduce)/2),i);
    imagesc(to,fo,vec);
    axis xy;
    set(gca,'xticklabel','','yticklabel','')
    if i==1
        title(['Every ' num2str(reduce) ' calls'])    
    end
end