#include <stdio.h>
#include <stdlib.h>


int i;
short int buffer[MAXSIZE];
FILE *in;
FILE *out;
char ifile[20],ofile[20];
/**/
main(){
for (i = 0; i < 5000; i++)
{
	sprintf(ifile,"C:\\Document\\file%04d.dat",i);
/*Since i goes from 0-5000, file0000.dat to file5000.dat will be written into ifile*/
	sprintf(ofile,"C:\\Document\\o_file%04d.dat",i);
/*Similar to ifile, o_acq0000.dat to acq0010.dat will be written to ofile*/
	
	in = fopen(ifile,"rb");
	out = fopen(ofile,"wb");
	fread(buffer,sizeof(short int),MAXSIZE,in);
	fwrite(buffer,sizeof(short int),MAXSIZE,out);
	
	fclose(in);
	fclose(out);
	return(0);
}
}

