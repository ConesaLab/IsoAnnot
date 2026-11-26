***********************************************************************************
* NucImport FOR DISTRIBUTION						          *
* Copyright Ahmed Mehdi and Mikael Boden 2011					  *
*										  *	
* LICENCE 									  *
*										  *
* This program is free softwareware: you can redistribute it under the 		  *
* terms of the GNU General Public License as published by the			  *
* Free Software Foundation.					  *
*										  *	
* We hope that this program will be useful but WITHOUT ANY WARRANTY.		  *	
*										  *
* You should have received a copy of the GNU General Public License		  *
* along with this program (http://www.gnu.org/licenses/). 			  *	
*										  *	
***********************************************************************************

README file

In this text file we provide the basic installation steps and usage instructions.

INSTALLATION

We expect that the users are familiar with basic java commands (example is also given below) to run this software.
The user should place the all files under NucImport folder within a accessible folder for java commands.
The user should place their fasta files in the same folder.  

SYNTAX FOR USAGE

java (maximum memory heap size) -jar (name of jar file) (name of your fasta file) (training model) (name of species) (indication of ID)

Parameter 1: (maximum memory heap size)
-Java has a couple of settings that help control how much memory it uses:
  - Xmx sets the maximum memory heap size
  - Xms sets the minimum memory heap size
Our program requires a maximum heap size of 2GB

Parameter 2: (name of jar file)
- The name of our jar file is NucImport.jar

Parameter 3: (name of your fasta file)
- Here you give your fasta file name. 

Parameter 4: (training model)
- We used six fold cross validation. Therefore we have provided six training models. Provide the name of training model (e.g. Model1 or Model2...or Model6)

Parameter 5: (name of species)
- We have tested our model on two species, YEAST and MOUSE. Provide the name of species by indicating Yeast or Mouse.

Parameter 6: (indication of ID)
- Our model provides optimal results if Ensemble IDs are provided for mouse and ORF IDs are provided for yeast. Please indicate by writing ID=T or ID=F
ID=T, means Ensemble mouse IDs are provided if mouse data is tested or ORF IDs are provided if yeast data is tested.
ID=F, means other IDs are provided. Our model considers the importins and Ran nodes to be unobserved if ID=F.

Examples:

1) java -Xmx2048m -jar NucImport.jar test.fa model6 Mouse ID=T

In case you already have more than 2G memory space, do not use Xmx2048m parameter (in such case an example is shown below),

2) java -jar NucImport.jar test.fa model6 Mouse ID=T


OUTPUT

The output of NucImport contains the protein ID followed by its import probability followed by its localization position and class (if any).



