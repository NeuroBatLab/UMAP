function [DATAindx]=get_random_calls(Pos,IDs,indxid,N,rep,weighted)
%grabs random # (=N) of calls in IDs
%this is done for certain number of repetitions (=rep)
%INPUT:
%       1) vector which defines how many rows (=snippets) one vocalization
%       had; given from prep_audio_umap.m
%       2) vector of IDs, same length as Pos; given from prep_audio_umap.m
%       3) individual bat IDs from who vocalizations should be grabbed
%       4) number of vocalizations to be grabbed
%       5) repetition number, how often random N calls should be grabbed
%       for each indxid
%       6) weighted (=1) indicates whether same number of calls should be
%       grabbed per individual or if 0, random number of calls independent
%       of individuals
%OUTPUT: indices for each repetition for the rows in the umap prepped matrix
%from prep_audio_umap.m
DATAindx=cell(1,rep);
for rp=1:rep
    indices=[];
    if weighted
        %grab N number of random calls per individual
        for bn=1:length(indxid)
            indx=find(IDs==indxid(bn));
            seq=randperm(length(indx));
            seq=seq(1:N);
            for i=1:N
                if indx(seq(i))-1>0
                    indices=[indices sum(Pos(1:indx(seq(i))-1))+1:...
                        sum(Pos(1:indx(seq(i))-1))+Pos(indx(seq(i)))];
                else
                    indices=[indices 1:Pos(indx(seq(i)))];
                end
            end
        end
    else
        %grab N number of calls independent of individual
        seq=randperm(length(Pos));
        seq=seq(1:N);
        for i=1:N
            if (seq(i))-1>0
                indices=[indices sum(Pos(1:(seq(i))-1))+1:...
                    sum(Pos(1:(seq(i))-1))+Pos((seq(i)))];
            else
                indices=[indices 1:Pos((seq(i)))];
            end
        end
    end
    DATAindx{rp}=indices;
end