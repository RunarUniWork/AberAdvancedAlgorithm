#include <stdio.h> /* needed for input, output, ... */
#include <stdlib.h> /* needed for EXIT_SUCCESS, ... */
#include <ctype.h> /* needed for isdigit() */
#include <string.h> /* needed for memset() */
#include <glpk.h> /* the linear programming toolkit */

#define SOURCENODE 0
#define SINKNODE 1

/* global variables -- YOU SHOULD NOT CHANGE THIS! */
/* (You are allowed to add your own if you want.) */
int numRows, numCols; /* number of rows and columns of the puzzle */
int *input; /* array holding the input */
int *solution; /* array holding the solution */
int debug; /* flag for debug mode; 1 means debug mode, 0 means debug off */

/* prototypes of functions -- YOU SHOULD NOT CHANGE THIS! */
/* (Feel free to add your own as you like.) */
int readInput(char *filename); /* reads puzzle from file */
/* readInput creates and fills the global variables as needed */
/* it returns 0 if all is okay and 1 otherwise */
int computeSolution(void); /* computes a solution if possible */
/* the return value is 1 if there is a solution and 0 otherwise */

void printPuzzle(int *puzzle);


/* This is the function that actually solves the problem. */
/* It is currently basically empty and not functional. */
/* Your own implementation needs to go in here. */
int computeSolution(void) {
	int totNodes = numRows * numCols + 2; //sum of nodes in grid + source & sink

	//======INPUT=2=GRAPH=======
	double upperbound = 1; //Sets the common upperbound for all edges

	int totEdges=0;
	int totDye=10;  
	int capacity =0; //Is for counting numbers of colours
	//Sett up graph
	double **graph; //Distance graph
	graph = malloc(sizeof(double *)*totNodes);
	//Allocate empty graph
	for(int row = 0; row < totNodes; row++){
		graph[row] = malloc(sizeof(double *)*totNodes);
		for(int col=0; col < totNodes; col++){ graph[row][col] = 0.0; }
	}

	//Fill graph with edge values
	for(int row = 0; row < numRows; row++)
		for(int col=0; col < numCols; col++){
			int node = row*numCols + col + 2;//convert [row][col] to node(2 first s&t)
			//Add outgoint edges that is connects node horizontally and vetically
			if (col != 0) 		{graph[node][node-1] = upperbound;		totEdges++;}
			if (col < numCols-1){graph[node][node+1] = upperbound;		totEdges++;}
			if (row != 0)		{graph[node][node-numCols] = upperbound;totEdges++;}
			if (row < numRows-1){graph[node][node+numCols] = upperbound;totEdges++;}
		}

	//Add source & sink nodes to graph
	int nodes[9] = { 0 }; //Array of max 9 values
   	int nodesFound = 0; //number of nodes that has been found
	//Go through the board to set >0 to source or sink
	for(int row = 0; row < numRows; row++)
		for(int col=0; col < numCols; col++)
			if(input[row*numCols+col] > 0){
				int found = 0; //Store if value is found
				for(int i=0; i < nodesFound; i++) //Check value already found?
					if(nodes[i] == input[row*numCols+col]){//Already found!
						found = 1;
						//Remove all edges going out of this node, except to the sink 
						for(int del=0;del<totNodes;del++)
							if(graph[row*numCols+col+2][del] != 0){
								graph[row*numCols+col+2][del] = 0;
								totEdges--;
							}
						graph[row*numCols+col+2][SINKNODE]=input[row*numCols+col];//Add path to sink
						totEdges++;
					}
				if(found == 0){//If not previously observed add edge from source
					nodes[nodesFound] = input[row*numCols+col];
					nodesFound++;
					//Remove all incomming edges to this edge, except from source
					for (int del=0; del < totNodes; del++)
						if(graph[del][row*numCols+col+2] != 0){
							graph[del][row*numCols+col+2]=0;
							totEdges--;
						}
					graph[SOURCENODE][row*numCols+col+2]=input[row*numCols+col];//Add path from source
					totEdges++;
				}
			}
	
	//======DO=LINEAR=PROGRAM============
	int edge = 0; //Counter for what edge is being observed

	glp_prob *lp; 
	lp = glp_create_prob();

	//Add edges to the lp problem with the bounds from graph
	glp_add_cols(lp, totEdges*totNodes*totDye);
	for(int row=0; row < totNodes; row++)
		for(int col=0; col < totNodes; col++)
			if(graph[row][col] > 0.0){
				for(int dye=1;dye<totDye;dye++){
					edge++;
					//If Source or sink, it only has cap for given colour
					if(row==SOURCENODE||row==SINKNODE||col==SOURCENODE||col==SINKNODE)
						if(graph[row][col] != dye) continue;
					glp_set_col_bnds(lp, edge, GLP_DB, 0.0, graph[row][col]);
				}
			}

	//Objective function(Maximize all edges out from source)
	glp_set_obj_dir(lp, GLP_MAX);
	edge=0;	
	for(int row=0; row < totNodes; row++)
		for(int col=0; col < totNodes; col++)
			if(graph[row][col] > 0.0)
				for(int dye=1; dye<totDye; dye++){
					edge++;
					if(row != SOURCENODE) continue;//Source is only to create flow
					if((int)graph[row][col] != dye) continue;
					glp_set_obj_coef(lp, edge, 1.0);
				}

	//Define constraints
	int constraintNumb = 1;
	glp_add_rows(lp, totNodes*2*totDye);//Set number of constraints for all transit nodes
	for(int node=0; node<totNodes; node++){
		//Sink and Source don't follow flow conservation
		//Set ^ into arrays for all diffconstraints ([2 +dye])
		//[1] Whole Flow conservation
		//[2] Sum of all out is <=1 
		//[3++] flow conservation for each dye
		int x = 2+totDye; //Short var for amount of diff constraint for each node
		int ruleIndex[x];
		int index[x][totNodes*totDye];
		double value[x][totNodes*totDye];
		for(int i=0; i<x; i++)
			ruleIndex[i] = 0;	

		for(int row = 0, edges=0; row < totNodes; row++)
			for(int col=0; col < totNodes; col++)
				if(graph[row][col] > 0.0)
					for(int dye=1; dye<totDye; dye++){
						edges++;
						if(node == SOURCENODE || node == SINKNODE) continue;
						if(col==node){ //Is it incomming node...
							ruleIndex[0]++;//Flow conservation
							index[0][ruleIndex[0]]=edges;
							value[0][ruleIndex[0]]=1.0;
	
							ruleIndex[dye+2]++; //Flow conservation for each dye
							index[dye+2][ruleIndex[dye+2]]=edges;
							value[dye+2][ruleIndex[dye+2]]=1.0;
						}else if(row==node){//...is it outgoing node
							ruleIndex[0]++;	//Flow conservation
							index[0][ruleIndex[0]]=edges;
							value[0][ruleIndex[0]]=-1.0;
						
							ruleIndex[1]++; //Sum of all out is >=0
							index[1][ruleIndex[1]]=edges;
							value[1][ruleIndex[1]]=1.0;

							ruleIndex[dye+2]++; //Flow conserv for each dye
							index[dye+2][ruleIndex[dye+2]]=edges;
							value[dye+2][ruleIndex[dye+2]]=-1.0;
						}
					}

		glp_set_row_bnds(lp, constraintNumb, GLP_FX, 0.0, 0.0); //Flow conservation
		glp_set_mat_row(lp, constraintNumb, ruleIndex[0], index[0], value[0]);
		constraintNumb++;
		glp_set_row_bnds(lp, constraintNumb, GLP_UP, 1.0, 1.0); //Sum of outgoing flow
		glp_set_mat_row(lp, constraintNumb, ruleIndex[1], index[1], value[1]);
		constraintNumb++;
		for(int dye=3; dye<totDye+2; dye++){ //Flow conservation 
			glp_set_row_bnds(lp, constraintNumb, GLP_FX, 0.0, 0.0);
			glp_set_mat_row(lp, constraintNumb, ruleIndex[dye], index[dye], value[dye]);
			constraintNumb++;
		}
	}
	
	//=======SOLVE=LP========
	glp_term_out(0);
	glp_simplex(lp, NULL);
	
	//Print out paths 
	if(debug == 1){ //Debugging method for printing the raw path
		printf("Max Flow is: %2.1f\n", glp_get_obj_val(lp));
		for(int row = 0, edges=0; row < totNodes; row++)
			for(int col=0; col < totNodes; col++)
				if(graph[row][col]>0.0){
					for(int dye=1; dye<totDye; dye++){
						edges++;
						double flow = glp_get_col_prim(lp, edges);
						if(flow > 0.0){
							if(row == SOURCENODE)
								printf("Source");
							else 
								printf("%d[%d][%d]",row-2,(row-2)/numCols,(row-2)%numCols);
							printf("\t-%d>\t", dye);
							if(col == SINKNODE)
								printf("Sink");
							else 
								printf("%d[%d][%d]",col-2,(col-2)/numCols,(col-2)%numCols);
							printf(" has flow: %f\n",flow);
						}
					}
				}
		printf("\n------\n");
	}
	//=======FLOW=PATH=2=MATRIX========
	for(int row = 0, edges=0; row < totNodes; row++)
		for(int col=0; col < totNodes; col++)
			if(graph[row][col]>0.0)
				for(int dye=1; dye<totDye; dye++){
					edges++;
					double flow = glp_get_col_prim(lp, edges);
					if(flow > 0.001){
						if(col <2) continue;
						int node = col-2;
						if(solution[node] !=0 && solution[node] != dye)
							capacity = 99; //Trying to cross lines, set check to big
						solution[node] = dye; //Set node to edges represented colour
					}
				}

	//Count amount of colours 
	for(int loop = SOURCENODE; loop < nodesFound; loop++){
	   if(debug == 1) printf("%d ", nodes[loop]); //Print colours
	   capacity++;
	}
	if(debug == 1){ //DEBUG to force print of input and found solution 
		printf(" <- is all colours found\n");
		printPuzzle(input);
		printf("\n");
		printPuzzle(solution);
		printf("\n");
	}

	//Cleanup
	for(int row = 0; row <totNodes; row++){free(graph[row]);}
	free(graph);

	if((int)glp_get_obj_val(lp) == capacity) return 1;
	return 0; /* this is not true for every puzzle, of course */
}

/* YOU SHOULD NOT CHANGE ANYTHING BELOW THIS LINE! */

/* printPuzzle(int *puzzle) prints either an input or a solution */
void printPuzzle(int *puzzle) {
    int i, j; /* loop variables to go over rows and columns */
    for ( i=0; i<numRows; i++ ) { /* go over rows */
      for ( j=0; j<numCols; j++ ) { /* go over columns */
        fputc('0'+puzzle[i*numCols+j], stdout); /* print the next char */
      }
      fprintf(stdout, "\n"); /* end the current line, start new one */
    }
}

int main(int argc, char **argv) {
	int i; /* used to run over the command line parameters */
	if ( argc<2 ) { /* no command line parameter given */
		fprintf(stderr, "Usage: %s [file1] [file2] [file3] [...]\n"
      "Where each [file] is the name of a file with a puzzle.\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	if ( argv[1][0]=='-' && argv[1][1]=='d' && argv[1][2]==0 ) {
    /* If the first parameter is -d we activate debug mode. */
		debug=1; /* switch debug mode on */
		fprintf(stdout, "DEBUG: Debug mode activated\n"); /* be explicit about it */
	} else {
		debug=0; /* switch debug mode off */
	}

  for ( i=1+debug; i<argc; i++ ) { /* go over remaining command line parameters */
    if ( readInput(argv[i]) ) { /* try to read file */
      /* returned with error message */
      fprintf(stderr, "%s: Cannot read puzzle with filename %s. Skipping it.\n",
        argv[0], argv[i]);
    } else { /* input read successfully */
			fprintf(stdout, "%s: Looking at the following puzzle:\n", argv[i]);
      if ( computeSolution() ) { /* compute a solution if one exists */
        fprintf(stdout, "%s: Found the following solution:\n", argv[i]);
        printPuzzle(solution);
      } else {
        fprintf(stdout, "%s: Puzzle has no solution\n", argv[i]);
      }
      /* free memory for next input */
      free(input);
      free(solution);
    }
  }
	return EXIT_SUCCESS;
}

/* checkFile(FILE *fh) performs basic checks and sets numRows/numCols */
/* return value 1 indicates an error; otherwise 0 is returned */
int checkFile(FILE *fh) {
  char c;
  int rows, cols; /* used to determine number of rows and columns */
  int read; /* counts number of digits read in the current row */
  int firstRow; /* indicates if we are reading the very first row */

  firstRow=1; /* we start in the first row */
  rows=cols=read=0;
  while ( !feof(fh) ) {
    c = fgetc(fh); /* read the next char from the file */
    if ( isdigit(c) ) { /* normal character read */
      read++; /* count the digit we just read */
      if ( ( !firstRow) && ( read>cols) ) {
        if ( debug ) {
          fprintf(stdout, "DEBUG: Row %d too long (%d, %d).\n", rows+1, read, cols);
        }
        return 1; /* flag error because row is too long */
      }
    } else {
      if ( ( c=='\n' ) || ( c==(char)-1 ) ) { /* end of line read */
        if ( read>0 ) {
          rows++; /* count the completed row if it was not empty */
        }
        if ( firstRow ) { /* very first row read */
          cols=read; /* accept number of characters as number of columns */
          firstRow=0; /* not in the first row anymore after this */
          if ( debug ) {
            fprintf(stdout, "DEBUG: %d columns per row expected\n", cols);
          }
        } else {
          if ( ( read>0 ) && ( read!=cols ) ) { /* rows too short */
            if ( debug ) {
              fprintf(stdout, "DEBUG: Row %d too short.\n", rows+1);
            }
            return 1; /* flag error because row is too short */
          }
        }
        read=0; /* reset number of characters in current row */
      } else { /* illegal character found */
        if ( debug ) {
          fprintf(stdout, "DEBUG: Illegal character %c found.\n", c);
        }
        return 1; /* stop reading because of the error */
      }
    }
  }
  if ( read>0 ) {
    rows++; /* last row was not ended with newline */
  }
  /* use the determined size and prepare for reading */
  numRows = rows;
  numCols = cols;
  rewind(fh); /* reset to the beginning of the file to read the actual input */
  return 0; /* signal all went well */
}

/* readInput(*char filename) reads the input and stores it */
/* return value 1 indicates an error; otherwise 0 is returned */
int readInput(char *filename) {
  int i, j; /* loop variables to go over the columns and rows of the input */
  int check[10]; /* array to check colours come in pairs */
  int keepReading; /* used to skip over newline */
  char c; /* next char */
  FILE *fh;

  if ( ( fh = fopen(filename, "rt") ) == NULL ) {
    return 1;
  }

  /* perform basic checks and determine size of puzzle */
  if ( checkFile(fh) ) { /* there was a problem */
    fclose(fh);
    return 1; /* signal error */
  }
  if ( ( input = (int *)malloc(sizeof(int)*numRows*numCols) ) == NULL ) {
    if ( debug ) {
      fprintf(stdout, "DEBUG: Unable to allocate %ld bytes of memory.\n",
        sizeof(int)*numRows*numCols);
    }
    fclose(fh);
    return 1;
  }
  if ( ( solution = (int *)malloc(sizeof(int)*numRows*numCols) ) == NULL ) {
    if ( debug ) {
      fprintf(stdout, "DEBUG: Unable to allocate %ld bytes of memory.\n",
        sizeof(int)*numRows*numCols);
    }
    free(input);
    fclose(fh);
    return 1;
  }
	memset(solution, 0, sizeof(int)*numRows*numCols);/*initialise solution empty*/
  if ( debug ) {
    fprintf(stdout, "DEBUG: Will read %dx%d sized puzzle\n", numRows, numCols);
  }
  /* prepare to count different digits */
  for ( i=0; i<10; i++ ) {
    check[i]=0;
  }
  /* Size is given in numRows, numCols; now we read */

  for ( i=0; i<numRows; i++ ) { /* go over rows */
    for ( j=0; j<numCols; j++ ) { /* go over columns */
      do {
        keepReading=1; /* prepare to skip over newline */
        c = fgetc(fh); /* get next digit */
        if ( isdigit(c) ) { /* store and count digit */
          input[i*numCols+j]=(int)(c-'0'); /* convert char to int */
          check[input[i*numCols+j]]++;
          keepReading=0; /* mark digit as read */
        }
      } while ( keepReading );
    }
  }
  for ( i=1; i<10; i++ ) {
    if ( ( check[i]!=0 ) && ( check[i]!=2 ) ) {
      if ( debug ) {
        fprintf(stdout, "DEBUG: Colour %d appears %d times.\n", i, check[i]);
        printPuzzle(input);
      }
      free(input);
      free(solution);
      fclose(fh);
      return 1;
    }
  }

  fclose(fh); /* close file after reading the input */
  return 0; /* signal all went well */
}
