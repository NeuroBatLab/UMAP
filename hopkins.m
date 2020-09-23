function H = hopkins (dataIn,xlim1, xlim2, ylim1, ylim2, m)
% HOPKINS calculates H for Hopkins statistics estimnation
%H = hopkins (dataIn,xlim1, xlim2, ylim1, ylim2, m)
sizeData = size(dataIn);
nPoints = sizeData(1); %number of datapoints
sizeVec = [xlim2-xlim1, ylim2 - ylim1];
%generate m random sampling points (not data points)
% sp(:,1) = xlim1 + sizeVec(1)*rand(m,1);
% sp(:,2) = ylim1 + sizeVec(2)*rand(m,1);
for i=1:2
sp(:,i)=normrnd(nanmean(dataIn(:,i)),nanstd(dataIn(:,i)),m,1);
end
%select m random data points:
pointIndex = [];
sdp = [];
ii = 1;
while ii <= m
    rn = round(nPoints * rand);
    if ~any(pointIndex == rn)
        if ~rn==0
            pointIndex(ii) = rn;
            sdp(ii,:) = dataIn (rn,:);
            ii = ii + 1;
        end
    end
end
% minimal distance between random sampling points and data points:
U = sum(min(dist2 (dataIn, sp), [], 1)); %minimum distances
% U=sum(min(pdist2(dataIn,sp),[],1));

leaveIndex = removerows((1:nPoints)',pointIndex); %excludes data points in sdp
dataR = dataIn(leaveIndex,:);
W = sum(min(dist2 (dataR, sdp), [], 1));

H = U/(U+W);