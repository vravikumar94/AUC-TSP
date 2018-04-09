
#include <stdio.h>
#include <string.h>
#include <stdlib.h>         
#include <strings.h>
#include  <math.h>
#include <stdbool.h>         

/////////////////////////
////	Ovarian dataset
/////////////////////////

/*
	#define MAX_GENES 1535
	#define Subject 53
	#define Sick 30
	#define Healthy 23
*/

/////////////////////////
////	Colon dataset
/////////////////////////

/*
	#define MAX_GENES 2000
	#define Subject 62
	#define Sick 22
	#define Healthy 40
*/

/////////////////////////
////	Leukemia dataset
/////////////////////////

#define MAX_GENES 3571
#define Subject 72
#define Sick 47
#define Healthy 25


/////////////////////////
////	Breast-LN dataset
/////////////////////////

/*
#define MAX_GENES 7129
#define Subject 49
#define Sick 25
#define Healthy 24
*/

/////////////////////////
////	Breast-ER dataset
/////////////////////////

/*
#define MAX_GENES 7129
#define Subject 49
#define Sick 24
#define Healthy 25
*/

/////////////////////////
////	DLBCL dataset
/////////////////////////

/*
#define MAX_GENES 7129
#define Subject 58	
#define Sick 32
#define Healthy 26
*/

/////////////////////////
////	DLBCL-FL dataset
/////////////////////////

/*
#define MAX_GENES 7129
#define Subject 77
#define Sick 58
#define Healthy 19
*/


/////////////////////////
////	Prostate dataset
/////////////////////////

/*
#define MAX_GENES 12533
#define Subject 102
#define Sick 52
#define Healthy 50
*/


float highdisease[MAX_GENES][MAX_GENES]={0}, lessdisease[MAX_GENES][MAX_GENES]={0}, highnormal[MAX_GENES][MAX_GENES]={0}, lessnormal[MAX_GENES][MAX_GENES]={0};
float Delta[MAX_GENES][MAX_GENES];
float data[MAX_GENES][Subject];
float a;
int NumberOfSubjects;
int label[Subject];
char* GeneID[MAX_GENES];
float ss,sd,sf,sa;
float Max[MAX_GENES][MAX_GENES] = {0};
float pp1[MAX_GENES]={0},pp2[MAX_GENES]={0};
float FirstMax = 0;
float FirstMaxAUC = 0;

/////////////////////////
////	Read the data file
/////////////////////////

void Read(FILE *file){

char c[100]={0};
int i = 0,j;

while ( fscanf(file,"%s", c) != EOF && i < MAX_GENES)
{

if (strstr(c,"Hsa.")!=NULL||strstr(c,"_")!=NULL||strstr(c,"g")!=NULL)
{
if( i > MAX_GENES )
{printf("\ni =  %d\n",i);}

i++;
GeneID[i]=strdup(c);
if( i > MAX_GENES )
{printf("\ni =  %d\n",i);}
if( i > MAX_GENES )
{printf("index\n");}
j=0;
}

else
{
data[i][j]=atof(strdup(c));
j++;
}
}
}

/////////////////////////
////	label the classes
/////////////////////////

void Label(){
int i;
for (i=0; i<Subject; i++)
{
if (i<Healthy)
{
label[i]=0;
}
else
label[i]=1;
}
}

/////////////////////////
////	Calculate the Delta_i_j for all the possible Gene pairs 'i' and 'j' 
/////////////////////////

float delta(float sa, float ss, float sd, float sf){
int i,j,k,l;
int s;
float p1,p2;

p1=sa/(sa+sd);
p2=ss/(ss+sf);
a= fabs(p1-p2);   

return a;
}

/////////////////////////
////	Find the maximum Delta value (the TSP score) 
/////////////////////////

float Maximum(){
int k, i, j;
float temp;
float max[MAX_GENES];
int a;

for (k=0; k<MAX_GENES; k++){
max[k] = 0;
for (i = 0; i < (MAX_GENES); i++){
if (Delta[k][i] > max[k]){
max[k] = Delta[k][i];
}
}
}

temp=0;
for (i = 0; i < (MAX_GENES); i++){
if ( max[i] > temp){
temp=max[i];
}
}

return temp;
}

/////////////////////////
////	Count number of times Gene 'i' < Gene 'j' in the healthy  and 
////						  Gene 'i' < Gene 'j' in the diseased groups 
////				Using the TSP method
/////////////////////////

void Count_TSP(int NumberOfSubjects){

int i,j,s,k;

for (i=0; i<MAX_GENES; i++){			
for (j=0; j<MAX_GENES; j++){
for (s=0; s<NumberOfSubjects; s++){
if(data[i][s]<data[j][s] && label[s]==1){
lessdisease[i][j]++;
}
if(data[i][s]<data[j][s] && label[s]==0){
lessnormal[i][j]++;
}
if(data[i][s]>=data[j][s] && label[s]==1){
highdisease[i][j]++;
}
if(data[i][s]>=data[j][s] && label[s]==0){
highnormal[i][j]++;
}
}
}
}

for (i=0; i<MAX_GENES; i++){						
for (j=0; j<MAX_GENES; j++){
Delta[i][j]= delta(lessdisease[i][j],lessnormal[i][j], highdisease[i][j], highnormal[i][j]);
}
}

}

/////////////////////////
////	Count number of times Gene 'i' < Gene 'j' in the healthy  and 
////						  Gene 'i' < Gene 'j' in the diseased groups 
////				Using the AUCTSP method
/////////////////////////

void Count_AUCTSP(int NumberOfSubjects){

int i,j,s,ss,k;

for (i=0; i<MAX_GENES; i++){			
for (j=i+1; j<MAX_GENES; j++){
for (s=0; s<Healthy; s++){
for (ss=0 ; ss<Healthy; ss++){
if(data[i][s]<data[j][ss])
lessnormal[i][j]++;
else if(data[i][s]>=data[j][ss])
highnormal[i][j]++;
}
}
}
}

for (i=0; i<MAX_GENES; i++){			
for (j=i+1; j<MAX_GENES; j++){
for (s=Healthy; s<Subject; s++){
for (ss=Healthy; ss<Subject; ss++){
if(data[i][s]>=data[j][ss])
highdisease[i][j]++;
else if(data[i][s]<data[j][ss])
lessdisease[i][j]++;
}
} 
}
}


for (i=0; i<MAX_GENES; i++){								
for (j=i+1; j<MAX_GENES; j++){
Delta[i][j]= delta(lessdisease[i][j],lessnormal[i][j], highdisease[i][j], highnormal[i][j]);
}
}
}

