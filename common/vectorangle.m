function angle = vectorangle(a,b)

angle = acos(dot(a,b,1)./(norm(a)*norm(b)));

end