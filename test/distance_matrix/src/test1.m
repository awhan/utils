load d.out;
load c.out;
load z.out;
d=d(:,2:end);
c=c(:,2:end);
z=z(:,2:end);
dm = sqrt(squeeze(sumsq(d - permute(c, [3 2 1]), 2)));
norm(dm-z)
