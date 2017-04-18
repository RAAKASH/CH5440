function [Fu] = Obj(B,P,Z)
m = 7;
n=3;
Fu = rms(rms((Z-(reshape(B,m,n)*P)))).^2;

end

