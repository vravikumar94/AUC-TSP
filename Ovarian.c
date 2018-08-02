#include "Functions.h"



int main (int argc, char *argv[]) { 
	
	
	
	int i,j,s; 
	int k,l;
	float test,train;
	int Sign = 0;
	FILE *TSP = fopen("TSP.txt", "w");
	FILE *AUCTSP = fopen("AUCTSP.txt", "w");
	
	FILE *file;
	
	if(argc!=2)
  	{
     	printf("Insufficient arguments\n");
     	exit(0);
	}
	file = fopen(argv[1],"r");
	
	
	if (AUCTSP == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}
	
	if (TSP == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}
	//	Read the data file
	Read(file);
//	Label the healthy and diseased classes
Label();
//	Count number of times Gene 'i' < Gene 'j' in the healthy  and 
////			        Gene 'i' < Gene 'j' in the diseased groups 
////				Using the AUCTSP method 
Count_AUCTSP(Subject);
//	Find the Top Score 	
FirstMax = Maximum();
//	Print out the information of the Top Scoring Pair using the AUCTSP	
for (i=0; i<MAX_GENES; i++){
for (j=0; j<MAX_GENES; j++){
if (Delta[i][j] == FirstMax){
fprintf(AUCTSP,"Score\t\tGene 1\t\tGene 2\n");
fprintf(AUCTSP,"%f\t%s\t%s\n",Delta[i][j], GeneID[i], GeneID[j]);
fprintf(AUCTSP,"P[WT] = %f\t P[Mutant] = %f\n", lessdisease[i][j]/(lessdisease[i][j] + highdisease[i][j]),lessnormal[i][j]/(lessnormal[i][j] + highnormal[i][j]));


}

}

}

//	Count number of times Gene 'i' < Gene 'j' in the healthy  and 
////			        Gene 'i' < Gene 'j' in the diseased groups 
////				Using the TSP method 
Count_TSP(Subject);
//	Find the Top Score	
FirstMax = Maximum();
//	Print out the information of the Top Scoring Pair using the TSP
for (i=0; i<MAX_GENES; i++){
for (j=0; j<MAX_GENES; j++){
if (Delta[i][j] == FirstMax){
fprintf(TSP,"Score\t\tGene 1\t\tGene 2\n");
fprintf(TSP,"%f\t%s\t%s\n",Delta[i][j], GeneID[i], GeneID[j]);
fprintf(TSP,"P[Sick] = %f\t P[Healthy] = %f\n", lessdisease[i][j]/(lessdisease[i][j] + highdisease[i][j]),lessnormal[i][j]/(lessnormal[i][j] + highnormal[i][j]));


}

}

}

fclose(TSP);
fclose(AUCTSP);

} // end of the main method
