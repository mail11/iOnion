function g = SUMSQIntTL(Dk)

load TripleLayerProfile.txt;
data = TripleLayerProfile;
g =  sum((data-CNIntTL(Dk)).^2);


end
    