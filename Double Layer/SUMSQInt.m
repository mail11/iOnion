function g = SUMSQInt(Dk)

load DBLayerProfile.txt;
data = DBLayerProfile;
g =  sum((data-CNInt(Dk)).^2);


end
    