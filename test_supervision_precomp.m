clear all;
close all;
Train=rand(3e3,3e3);
for i=1:size(Train,1)
    Train(i,i)=0;
    for ii=1:size(Train,1)
        Train(ii,i)=Train(i,ii);
    end
end
Train(:,end+1)=0;
Train(1.8e3:end,end)=1;

Test=rand(3e3,3e3);
for i=1:size(Test,1)
    Test(i,i)=0;
    for ii=1:size(Test,1)
        Test(ii,i)=Test(i,ii);
    end
end
labs=[];
for f=1:size(Train,2)-1
    l=['f' num2str(f)];
    labs{f}=l;
end
%run template=train data-set
[reduction, umap,clustering]=run_umap(Train,'parameter_names',labs,...
    'label_column',size(Train,2),'metric','precomputed',...
    'verbose','graphic');
save([pwd '\umap_data\UMAP_saldeaftrain_dtw.mat'], 'umap','reduction','clustering','-v7.3');

[reductiont, umapt,clustering]= run_umap(Test, 'template_file',...
    [pwd '\umap_data\UMAP_saldeaftrain_dtw.mat'],'parameter_names',labs,...
    'metric','precomputed','verbose','graphic');

