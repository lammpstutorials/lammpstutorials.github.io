file=fopen('metadata.dat','wt');
for a=1:24
X=['./position.',num2str(a),'.dat ',num2str(-20+0.5*(a-1)),' 0.5'];
fprintf(file,X);
fprintf(file,'\n');
end
