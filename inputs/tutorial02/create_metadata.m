file=fopen('metadata.dat','wt');
for a=1:67
	X=['./position.',num2str(a),'.dat ',num2str(a*3-102),' 0.0205'];
	fprintf(file,X);
	fprintf(file,'\n');
end
