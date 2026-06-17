function weixing2(XX,dire,xuanzhuan)
alpha = 0.9;
d = 30;
yanse = [220 220 220];
yanse = yanse/norm(yanse);
mian1x = [XX(1)-d XX(1)+d XX(1)+d XX(1)-d];
mian1y = [XX(2)+d XX(2)+d XX(2)-d XX(2)-d];
mian1z = ones(1,4)*(XX(3)+d);
s = patch(mian1x,mian1y,mian1z,yanse,'FaceAlpha',alpha);
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

mian2x = [XX(1)-d XX(1)+d XX(1)+d XX(1)-d];
mian2y = [XX(2)+d XX(2)+d XX(2)-d XX(2)-d];
mian2z = ones(1,4)*(XX(3)-d);
s=patch(mian2x,mian2y,mian2z,yanse,'FaceAlpha',alpha);
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

mian3x = [XX(1)-d XX(1)+d XX(1)+d XX(1)-d];
mian3z = [XX(3)+d XX(3)+d XX(3)-d XX(3)-d];
mian3y = ones(1,4)*(XX(2)+d);
s=patch(mian3x,mian3y,mian3z,yanse,'FaceAlpha',alpha);
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

mian4x = [XX(1)-d XX(1)+d XX(1)+d XX(1)-d];
mian4z = [XX(3)+d XX(3)+d XX(3)-d XX(3)-d];
mian4y = ones(1,4)*(XX(2)-d);
s=patch(mian4x,mian4y,mian4z,yanse,'FaceAlpha',alpha);
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

mian5y = [XX(2)-d XX(2)+d XX(2)+d XX(2)-d];
mian5z = [XX(3)+d XX(3)+d XX(3)-d XX(3)-d];
mian5x = ones(1,4)*(XX(1)+d);
s=patch(mian5x,mian5y,mian5z,yanse,'FaceAlpha',alpha);
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

mian6y = [XX(2)-d XX(2)+d XX(2)+d XX(2)-d];
mian6z = [XX(3)+d XX(3)+d XX(3)-d XX(3)-d];
mian6x = ones(1,4)*(XX(1)-d);
s=patch(mian6x,mian6y,mian6z,yanse,'FaceAlpha',alpha);
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

yi1x = [XX(1)+1.5*d XX(1)+4*d XX(1)+4*d XX(1)+1.5*d];
yi1y = [XX(2)+0.5*d XX(2)+0.5*d XX(2)-0.5*d XX(2)-0.5*d];
yi1z = ones(1,4)*(XX(3));
s=patch(yi1x,yi1y,yi1z,yanse,'FaceAlpha',1);
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

yi2x = [XX(1)-1.5*d XX(1)-4*d XX(1)-4*d XX(1)-1.5*d];
yi2y = [XX(2)+0.5*d XX(2)+0.5*d XX(2)-0.5*d XX(2)-0.5*d];
yi2z = ones(1,4)*(XX(3));
s=patch(yi2x,yi2y,yi2z,yanse,'FaceAlpha',1);
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

x1x = [XX(1)+d XX(1)+1.5*d];
x1y = [XX(2) XX(2)-0.5*d];
x1z = [XX(3),XX(3)];
s=line(x1x,x1y,x1z,'Color','black');
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

x2x = [XX(1)+d XX(1)+1.5*d];
x2y = [XX(2) XX(2)+0.5*d];
x2z = [XX(3),XX(3)];
s=line(x2x,x2y,x2z,'Color','black');
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

x3x = [XX(1)-d XX(1)-1.5*d];
x3y = [XX(2) XX(2)-0.5*d];
x3z = [XX(3),XX(3)];
s=line(x3x,x3y,x3z,'Color','black');
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

x4x = [XX(1)-d XX(1)-1.5*d];
x4y = [XX(2) XX(2)+0.5*d];
x4z = [XX(3),XX(3)];
s=line(x4x,x4y,x4z,'Color','black');
rotate(s,dire,xuanzhuan,[XX(1),XX(2),XX(3)]);
hold on

end