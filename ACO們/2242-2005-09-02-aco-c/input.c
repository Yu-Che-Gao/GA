
#include"dec.h"
#include"func.h"
	float midht;

void read_data(char * fname)
{

	FILE *fin = fopen(inpfname , "r");
	char finp[strlen(inpfname)];
	strcpy(finp, inpfname);
	FILE *fgraph = fopen( strcat(finp, ".graph.out" ), "w");
	if(fin==NULL) printf("error in opening input file\n");

	int vno , valueno ,maxnv=0 ;
	#if (taco ==1)
	printf("enter the tournament size \n"); 
	fscanf(fin,"%d", &tourneysize);
	printf("tournament size set to %d \n\n",tourneysize); 
	#endif
	for(vno=0;vno<nov;vno++)
	{
		printf("enter the name for variable no %d \n", vno+1);
		fscanf(fin,"%s",(var[vno]).name);
		printf("enter the lower limit for %s \n",var[vno].name);
		fscanf(fin,"%f",&(var[vno]).llimit);
		printf("enter the upper limit for %s \n",var[vno].name);
		fscanf(fin,"%f",&(var[vno]).ulimit);
		printf("enter the increment for %s \n",var[vno].name);
		fscanf(fin,"%f",&(var[vno]).increment);
		
		find_nvalue(vno );
		
		(var[vno]).trail = (double *) calloc((var[vno]).nvalues ,sizeof(double));
		(var[vno]).prob = (double *) calloc((var[vno]).nvalues ,sizeof(double));
		if((var[vno]).nvalues > maxnv)
		{
			maxnv = (var[vno]).nvalues;
		}
	
	}//end of nov loop
	
	//midht = (maxnv  )/ 2.0; 
	midht = 0.0;
	for(vno=0;vno<nov;vno++)
	{
		for(valueno=0;valueno< (var[vno]).nvalues ; valueno++)
		{
		fprintf(fgraph , "%15d",vno+1);
		fprintf(fgraph , "%15f\n", midht - valueno + ((var[vno]).nvalues - 1.0) / 2.0 );
		}
		fprintf(fgraph,"\n");
	}
		fprintf(fgraph , "%15d",1 + nov);
		fprintf(fgraph , "%15f\n\n", midht );
	printf("%d is the maxnv and midht is %f \n\n",maxnv, midht);
	
	fclose(fin);
	fclose(fgraph );



	printf("tournament size set to %d \n\n",tourneysize); 
	//scanf("%d", &tourneysize);

}//end of read data

void print_data()
{

char finp[strlen(inpfname)];
strcpy(finp, inpfname);
FILE *input= fopen( strcat(finp, ".var.tex" ), "w");
printf("in print data finp is %s\n",finp); 
printf("before var.out in final stats inpfname is %s\n",inpfname); 

printf("var.out after  in final stats finp is %s\n",finp); 
printf("var.out fater  in final stats inpfname is %s\n",inpfname); 
int vno;
fprintf(input,"\\begin{table}[!hp] \\label{varlimits}\n  \\begin{center}\n \\caption{Variable details  for %s }\n \\begin{tabular}{rlr@{--}lr} \\hline \n", inpfname);
fprintf(input,"$i$& $\\mathcal{X}_i$& \\multicolumn{2}{c}{$\\mathcal{D}_i$} & $\\Delta_i  $  \\\\ \\hline \n");
for(vno=0;vno<nov;vno++)
{
fprintf(input,"%d& %s &%2.0f&  %2.0f & %2.0f \\\\\n",vno+1, (var[vno]).name, ((var[vno]).llimit),(var[vno]).ulimit,(var[vno]).increment);
}
fprintf(input,"\\hline \n \\end{tabular}\n \\end{center}\n \\end{table}\n");


fclose(input);

}//end of print data


void find_nvalue(int vno)
{
float nv,quotient;
int divisor;

nv =  ( ((var[vno]).ulimit - (var[vno]).llimit )) / ((var[vno]).increment );

divisor = (int) nv;

(var[vno]).nvalues =1 + divisor ; 

quotient = nv /(( float)divisor);

if(quotient !=1  );
{
 	if((fabs((quotient - 1.000000000000000) ) > 0.01))
	{
 	((var[vno]).nvalues) ++;
	}

 }

}//end of find n value

