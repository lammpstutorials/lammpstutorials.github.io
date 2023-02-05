file=fopen('metadata.dat','wt');
for a=1:50
	X=['./position.',num2str(a),'.dat ',num2str(a-25),' 1.5'];
	fprintf(file,X);
	fprintf(file,'\n');
end
