function c = randchar(n)

u = rand(1,n);
u = floor(u*26)+65;
c=char(u);