function [newtrindx]=check_dataforumap(DATA,Pos,minN_umap,IDs,batid,randprop,N)
newtrindx{1}=[];
%if not enough data for umap
if length(DATA)<minN_umap
    disp('N too small for this rep! Redoing random data set generation independent of individual')
    callN=round(randprop*length(IDs));
    [newtrindx]=get_random_calls(Pos,IDs,batid,callN,1,0);
    %if not enough data for umap
    if length(newtrindx{1})<minN_umap
        disp('N too small! Grabbing more data')
        callN=round(.8*length(IDs));
        [newtrindx]=get_random_calls(Pos,IDs,batid,callN,1,0);
        %if not enough data for umap
        if length(newtrindx{1})<minN_umap
            disp('N too small! Grabbing all data')
            newtrindx{1}=1:N;
        end
    end
end