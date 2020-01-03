/* version for sorting out normal tissue consumption, will use small domain for computational efficiency 
have to ensure that normal cells (0) consume at basal rate and then sort what the hell that is
this version HAS NECROTIC CELLS, HAS NORMAL CELLS, MM kinetics on oxygen consumption, no update on vessels
new version here has HashTable trying to calculate phylogenies with mutation (neutral)

BIT STRINGS FOR MUTATION AND DEATH

note oxygen has been toggled off by commenting line 290
*/

import java.io.*;
import java.util.*;
import static java.lang.Math.*;

public class CA3 {
	public static String version="8th dec 2015";
   	public boolean simulationFinished=false;
	MersenneTwisterFast random;
	public int mag = 1;
	public float [][] Oxygen = null; 
	public float [][] consumption = null; // make matrix for consumption
	public int [][] Age;
	int [][] Cells = null;  // Cells  and the cell type
					// 0 -> Healthy Cell
					// 1 -> Stem Cells
					// 2 -> Progenitor Cells
					// 3 -> Mature/Differentiated cells
					// 4 -> Necrosis

	//int [][] Vasculature = null; // Initial vasculature as a boolean-like lattice
	Bag cellList = null; // Used in the iterateCells method. Defined here for performance issues
	boolean finishedRun = false;
	int mutationNum = 1;  // counter for mutations, starts with 1 which represents the mutation of the first cancer cell
	float mutfreq = 0.01f; // frequency of mutation
	float deathprob = 0.01f;
	StringBuilder [][] carriedGenome = null; // instantiating carriedGenome
	StringBuilder initGenome;
	int lengthGenome = 1000;
	//Hashtable<Integer, Integer> tree;
	//Hashtable<Integer, Integer> timeTree;
	StringBuilder genomeToMod = new StringBuilder();
	//int [] mutTimePair = null;  //******
	
        // Vascular initialization
		//int spacing = 12; // for equal vascular spacing
    	//int vesselNumber = 50; // number of vessels to seed in alternative formulation
		int size = 100; //8*spacing+1; // Size of the system lattice
        int centre = size/2; //find centre of lattice
		int timestep = 0; // Current Number of timesteps in the simulation
        //int o2perTS = 2304; // 230400 Times we iterate oxygen per cell time step
		int dataWriteStartTime = 1; // when we start writing down the data
        int dataWriteStep = 500; // how often we write down the data
        int dataReportStep = 500; // how often we print the data
        //float consumption_rate = 0.0f; //usually 1


	// Model Parameters
        //float kDe = 0.1f; // 0.001728f;
        float hypoxia = 0.1f; // Hypoxia threshold
		float prolifThreshold = 0.11f; // minimal oxygen to proliferate
        //float Km = 0.01f; // Michaelis-Menten kinetic parameter
		float initOxygen = 1.0f;
		float[] proliferation = {
		0.025f,  /* healthy cells */ 
		0.05f, /* stem cells*/
		0.05f, /* progenitor cells */
		0.0f, /*mature/differentiated cells */
		0.0f, /*necrotic area */ // added this cell type
						};
/*
        float[] consumptionBasal = 
        { // 0.000375 and 0.00075
		0.0005f,  /* healthy cells 
		0.001f, /* stem cells
		0.001f, /* progenitor cells 
		0.001f, /*mature/differentiated cells 
		0.0f, /*necrotic area  // added this cell type
		};
	

    	float[] consumptionDivision = 
    	{
		5 * consumptionBasal[0], // The Oxygen consumption for normal cells that proliferate
		5 * consumptionBasal[1], // The Oxygen consumption for TICs that proliferate
		5 * consumptionBasal[2], // The Oxygen consumption for TACs that proliferate
		5 * consumptionBasal[3], // The Oxygen consumption for TDs that proliferate
		5 * consumptionBasal[4], // The Oxygen consumption for necrotic that proliferate which is zero
    	};
*/

	int maxProDivisions;
	int maxMatureCellAge;
	float asymmetricRatio;
	float pMotility=0.00f; // COMPLETELY ARBITRARY INDEED, REVISIT
	//float densityVasculature=0.004f; // 1/250
	//float avgVesselRad=2;
    boolean radiotherapy=false;
	
	// Measurement and STATISTICS
	int births=0;
	int deaths=0;
    int[][] stemBirthCounter = null; // tracks stem cell mitotic age - resets to zero with cell death
	int[][] stemDeathCounter = null; // tracks total stem deaths at a position over time
	int[][] stemBirthsTotal = null; // tracks total stem births at a position over time
	//int[][] TACBirthCounter = null; // tracks total TAC births at a position over time
	//int[][] TACDeathCounter = null; // tracks total TAC deaths at a position over time
	int[][] carriedmutation = null; // tracks most recent parental mutation at position [i][j] over time

/*
	// ******
	public CA3 (float DiffusionCoeff, float ConsumptionRate, float MichMen, float Hypox)
		{
		int time = (int) System.currentTimeMillis();
		random = new MersenneTwisterFast (time);
		kDe = DiffusionCoeff;
		consumption_rate = ConsumptionRate;
		Km = MichMen;
		hypoxia = Hypox;
		for (int i=0; i<5; i++) {
  		consumptionBasal[i] = consumptionBasal[i] * consumption_rate;
			}
		prolifThreshold = hypoxia + 0.01f; // minimal oxygen to proliferate
		//asymmetricRatio=SCSymmetricDivBalance;
		//maxProDivisions=maxDivisions;
		maxMatureCellAge=1;
		//densityVasculature=densityV;
		reset();
		resetVasculature();
	}


	public class Pair<X, Y> { 
  		public final X x; 
  		public final Y y; 
  		public Pair(X x, Y y) { 
   		 this.x = x; 
   		 this.y = y; 
 		 } 
		} 
*/

	public CA3 (float SCSymmetricDivBalance, int maxDivisions, float densityV)
	{
		int time = (int) System.currentTimeMillis();
		random = new MersenneTwisterFast (5); //choose a seed a number or 'time'
		asymmetricRatio=SCSymmetricDivBalance;
		maxProDivisions=maxDivisions;
		maxMatureCellAge=1;
		//densityVasculature=densityV;
		carriedGenome = new StringBuilder [size][size];
		//tree = new Hashtable<Integer,Integer>();
		//timeTree = new Hashtable<Integer, Integer>();
		reset();
		//resetVasculature();
	}


	// ******
	public void reset ()
	{
		//consumption = new float [size][size];		
		//Oxygen = new float [size][size];
		Age = new int [size][size];
		stemBirthCounter = new int [size][size];
		stemBirthsTotal = new int [size][size];
		stemDeathCounter = new int [size][size];
		//TACDeathCounter = new int [size][size];
		//TACBirthCounter = new int [size][size];
		StringBuilder initGenome = new StringBuilder();
		for (int i=0;i<lengthGenome;i++) {initGenome.append("0");};
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++) {
				carriedGenome[i][j] = (StringBuilder) initGenome;
				//carriedGenome[i][j] = (StringBuilder) initGenome;
				//Oxygen[i][j]=initOxygen;
				Age[i][j]=0;
				stemBirthCounter[i][j]=0;
									}

		resetCells();
	}

    // *******initialize******* domain with vessels
/*
	void resetVasculature ()
	{
	/* this routine for seeding an AMOUNT of vessels 
	    int k = 0; 
	     	    int i = 0;
	    	    int j = 0;
	     	    Vasculature = new int [size][size];
	     	    while (k < (vesselNumber)) {
	     		i = (int) Math.floor(Math.random()*size);
	    		j = (int) Math.floor(Math.random()*size);
	     		if (Vasculature[i][j]==0) {
	     		    Vasculature[i][j]=1;
	     		    k++;
	     		}
	     	    }
	    Vasculature[centre][centre]=1;


	    Vasculature = new int [size][size];	
 //use this routine to load vasculature from file
/*			
		try {
		  	FileInputStream inF = new FileInputStream("./text/LastVasc12");
		  	BufferedReader br = new BufferedReader(new InputStreamReader(inF));
		 
		  	String str = br.readLine();
		  	for (int i=0;i<=size;i++){
			  	StringTokenizer strTok = new StringTokenizer(str);
			  	for(int j=0; j<=size;j++){
			        Vasculature[i][j] = Integer.parseInt(strTok.nextToken());
			    }
			    str=br.readLine();
		  	}
		  			
		  	inF.close();
		} catch (IOException e) {
		System.out.println("dipshit");
		}

 //end grab file


	    	 for (int i=0;i<size;i++)
	         for (int j=0;j<size;j++)

		       	// if ( (i==0) || (i == (size-1) ) ||  (j==0) || ( j== (size-1) ) ) Vasculature[i][j]=1; // making the boundary vascularized
				//	if (i==0) Vasculature[i][j]=1; // making the left boundary vascularized
	    	    // if ((random.nextFloat()<=densityVasculature) && (Vasculature[i][j]==0)) Vasculature[i][j]=1; // randomly seeding probability from command line

		       if (((i==0) || (i%(2*spacing)==0)) && ((j+spacing)%(2*spacing)==0)) Vasculature[i][j]=1; // for equal spacing
		       else if (((j==0) || (j%(2*spacing)==0)) && ((i+spacing)%(2*spacing)==0)) Vasculature[i][j]=1; // for equal spacing

			//}


		try {             // write down the vasculature to a file for resuse if desired *********
                File dir = new File ("./text");
                dir.mkdir ();
		FileWriter outFile4 = new FileWriter("./text/LastVasc");
		PrintWriter outVasc = new PrintWriter(outFile4);
                for (int a=0;a<size;a++) {
                    for (int b=0;b<size;b++){
			outVasc.print(Vasculature[a][b]+", ");
                    }
                    outVasc.println("");
                }
		outVasc.close();
		} catch (IOException e) {
                e.printStackTrace();
                System.exit(-1);
		}  // cut here ************


	}
	*/
    // *******initialize******* domain with cells
	void resetCells ()
	{
		// int radius=50; dont think this is used
		Cells = new int [size][size];
		carriedmutation = new int [size][size];	
		
		// Fill the domain with healthy cells 0's first.
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++) {
				Cells[i][j]=0;carriedmutation[i][j]=0;
			}
		
		Cells[centre][centre] = 1; // initialize with a stem cell if you like
		carriedmutation[centre][centre] = 1; // the first cancer cell carries the tag '1' to start
        
        genomeToMod = new StringBuilder (carriedGenome[centre][centre]);
        genomeToMod.setCharAt(mutationNum-1,'1'); //  daughter carries new carriedGenome
        carriedGenome[centre][centre] = genomeToMod;

		//genomeToMod = carriedGenome[centre][centre];
	}

    // to measure euclidean distance (not yet used)
	final int distance (int x, int y, int i, int j)
	{
		double dis=Math.sqrt((x-i)*(x-i)+(y-j)*(y-j));
		return (int)Math.round(dis);
	}

    // This method makes the boundaries no flux ************

	final int[] convertCoordinatesNoFlux(int x, int y)
       {
               if (x < 0) x = 1;
               else if (x > size - 1) x = (size - 2);
               if (y < 0) y = 1;
               else if (y > size - 1) y = (size - 2);
               int[] result = new int[2];
               result[0] = x; result[1] = y;
               return result;
       }

     // This method makes the boundaries periodic **********

	final int[] convertCoordinates(int x, int y)
       {
               if (x < 0) x = size - 1;
               else if (x > size - 1) x = 0;
               if (y < 0) y = size - 1;
               else if (y > size - 1) y = 0;
               int[] result = new int[2];
               result[0] = x; result[1] = y;
               return result;
       }

	
	public void nextTimeStep ()
	{
	    // if (timestep==1000)   //zeros out vasculature at t=1000
        // for (int i=0;i<size;i++)
		// for (int j=0;j<size;j++) Vasculature[i][j]=0;
		births = 0;
		deaths = 0;
//		for (int i=0;i<o2perTS;i++) iterateOxygen();  // oh dear
		iterateCells();
		radiotherapy = false;
		
		//NEW
		int totalCells = 0;
		int totalHealthy = 0;
		int totalStem = 0;
		int totalProgenitor = 0;
		int totalMature = 0;
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++)
				if (Cells[i][j]<4) {
					totalCells++;
					if (Cells[i][j]==0) totalHealthy++;
					else if (Cells[i][j]==1) totalStem++;
					else if (Cells[i][j]==2) totalProgenitor++;
					else if (Cells[i][j]==3) totalMature++;
					else System.err.println ("wrong cell type");
				}
			
		// Not so new	
		if (timestep==0) System.out.println (+size+", "+mutfreq+", "+lengthGenome+", "+asymmetricRatio); // print (parent,child) pair

			//System.out.println ("% Timestep\t Cells\t Stem Cells \t Progenitor\t Mature");
		if (timestep%dataReportStep==0) {		
		    //System.out.println(timestep+"\t"+totalCells+"\t"+totalStem+"\t"+totalProgenitor+"\t"+totalMature+"\t"+births+"\t"+deaths+"\t"+((float)births/deaths));
		    //System.err.println(asymmetricRatio+" "+maxProDivisions+" "+totalCells+" "+totalStem);}
		    //System.err.println(mutationNum+" "+totalStem);
		}
		timestep++;
        
        // Finally let's write down the data
        if ( (timestep%dataWriteStep==0) && (timestep>dataWriteStartTime) ) {
            try {
                File dir = new File ("./text");
                dir.mkdir ();

		// cell matrix
		FileWriter outFile1 = new FileWriter("./text/cells"+timestep);
        PrintWriter outCells = new PrintWriter(outFile1);

        // attempt at hashtable
        //FileWriter outFileHm = new FileWriter("./text/pairs"+timestep); // new
        //PrintWriter outTable = new PrintWriter(outFileHm);  // new

                // carriedGenome
        FileWriter outFilecG = new FileWriter("./text/genomes"+timestep); // new
        PrintWriter outcG = new PrintWriter(outFilecG);  // new

                // attempt at hashtable
        //FileWriter outFileHm2 = new FileWriter("./text/timepairs"+timestep); // new
        //PrintWriter outTable2 = new PrintWriter(outFileHm2);  // new
		
		// oxygen matrix
                //FileWriter outFile2 = new FileWriter("./text/oxygen"+timestep);
                //PrintWriter outO2 = new PrintWriter(outFile2);

		// stem age matrix
                // FileWriter outFile3 = new FileWriter("./text/stemBirthCounter"+timestep);
                // PrintWriter outSBC = new PrintWriter(outFile3);

        // carried mutation matrix
        FileWriter outFile2 = new FileWriter("./text/carriedmutation"+timestep);
        PrintWriter outCM = new PrintWriter(outFile2);
/*
		// stem total birth matrix
                //FileWriter outFile4 = new FileWriter("./text/stemBirthsTotal"+timestep);
                //PrintWriter outSBM = new PrintWriter(outFile4);

		// stem death matrix
                //FileWriter outFile5 = new FileWriter("./text/stemDeathCounter"+timestep);
                //PrintWriter outSD = new PrintWriter(outFile5);

		// TAC birth matrix
                //FileWriter outFile6 = new FileWriter("./text/TACBirthCounter"+timestep);
                //PrintWriter outTB = new PrintWriter(outFile6);

		// TAC death matrix
                //FileWriter outFile7 = new FileWriter("./text/TACDeathCounter"+timestep);
                //PrintWriter outTD = new PrintWriter(outFile7);

        		//write hashtable
      			for (Integer key : tree.keySet()) {
      			outTable.print(+key+", "+tree.get(key)+"\r\n");
      			} //new
        		outTable.println(""); //new
        		outTable.close(); //new

        		//write hashtable
      			for (Integer value : tree.valueSet()) {
      			outTableNew.print(+tree.get(value)+", "+value+"\r\n");
      			} //new
        		outTableNew.println(""); //new
        		outTableNew.close(); //new


        		//write hashtable
        		for (Integer key : timeTree.keySet()) {
      			outTable2.print(+key+", "+timeTree.get(key)+"\r\n");
      			} //new
        		outTable2.println(""); //new
        		outTable2.close(); //new
*/
                for (int i=0;i<size;i++) {
                    for (int j=0;j<size;j++){
					outCells.print(Cells[i][j]+", ");
					outcG.print(carriedGenome[i][j]+", ");
                   		//outO2.print(Oxygen[i][j]+", ");
                        //outSBC.print(stemBirthCounter[i][j]+", ");
					outCM.print(carriedmutation[i][j]+", ");
						//outSBM.print(stemBirthsTotal[i][j]+", ");
						//outSD.print(stemDeathCounter[i][j]+", ");
						//outTB.print(TACBirthCounter[i][j]+", ");
						//outTD.print(TACDeathCounter[i][j]+", ");

                    }
                    	outCells.println("");
                    	outcG.println("");
                    	outCM.println("");
               	    	//outO2.println("");
                    	//outSBC.println("");
		    			//outSBM.println("");
						//outSD.println("");
						//outTB.println("");
						//outTD.println("");

                }
				outCells.close();
				outcG.close();
				outCM.close();
            	//outO2.close();
                //outSBC.close();
				//outSBM.close();
                //outSD.close();
                //outTB.close();
                //outTD.close();

            } catch (IOException e) {
                e.printStackTrace();
                System.exit(-1);
            }
        }
	}
	

    // main CELL CA loop************
	public boolean iterateCells()
 	{

/*
	// modify consumption matrix
	    	for (int i=0;i<size;i++)
	        for (int j=0;j<size;j++) consumption[i][j] = consumptionBasal[Cells[i][j]];
*/				
	//
		if (cellList==null) cellList = new Bag (size*size); 
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++) {
			    if (Cells[i][j] < 4 )  { // All tumour cell types have Cell > 0, now 0 corresponds to 'healthy cells' that consume at basal rate only 
					int[] p = new int[2];
					p[0] = i; p[1] = j;
					cellList.add(p);
				}
			}

		while (cellList.size() != 0) {
			// Select the next lattice element at random
			int randomElemIndex=0;
			if (cellList.size()>1) randomElemIndex = random.nextInt(cellList.size()-1);
	   		int[] point = (int[])cellList.get(randomElemIndex); 
		   	int rI = point[0];
		   	int rJ = point[1];

			cellList.remove(randomElemIndex); // Remove it from the cell list
		   	int cell = Cells[rI][rJ];
   
			// Cell death
			//if ((Oxygen[rI][rJ]<hypoxia)) {
			if ((random.nextFloat()<deathprob) && Cells[rI][rJ]>0){ // x% chances of dying 
				Age[rI][rJ]=0;
				if (Cells[rI][rJ]==1) stemDeathCounter[rI][rJ]++;
				if  (Cells[rI][rJ]<4 ) {
				//TACDeathCounter[rI][rJ]++;
				Cells[rI][rJ]=4; // was 0, now making necrotic area (truly empty)
				deaths++;
				stemBirthCounter[rI][rJ]=0;
				carriedmutation[rI][rJ] = 0; // empty space now has no mutations
				carriedGenome[rI][rJ] = initGenome; // empty space now has no mutations
				}}
			

			else if ((cell==3) && (Age[rI][rJ]>100*maxMatureCellAge)) { // added * to allow for an update each celltimestep/x **************************
				Age[rI][rJ]=0;
				Cells[rI][rJ]=4; // was 0, now making necrotic area (truly empty)
				//TACDeathCounter[rI][rJ]++;
				deaths++;
				}

			else if ( (radiotherapy) && (cell==2) && (random.nextFloat()>Oxygen[rI][rJ]) ) {
                		// Radiotherapy
                		Age[rI][rJ]=0;
				if (Cells[rI][rJ]==1) stemDeathCounter[rI][rJ]++;
				if ( (Cells[rI][rJ]==2) || (Cells[rI][rJ]==3) ) 
				//TACDeathCounter[rI][rJ]++;
				Cells[rI][rJ]=4; // make necrotic
				stemBirthCounter[rI][rJ]=0;
				deaths++;
            	}

			//healthy division
			else if ((cell==0) && (vacantSites(rI,rJ)>0)) {
			   if (proliferation[cell]>=random.nextFloat()) {// If tossing the coin we are to proliferate...
				    //if (Oxygen[rI][rJ]>prolifThreshold) { // AND the oxygen concentration is enough for division..
					//consumption[rI][rJ]=consumptionDivision[Cells[rI][rJ]];
				        int[] daughter = findEmptySite (rI,rJ);
						births++;
						Cells[daughter[0]][daughter[1]] = 0;
						carriedmutation[daughter[0]][daughter[1]] = 0; 
						carriedGenome[daughter[0]][daughter[1]] = initGenome;// resetting space to healthy cell with no mutations
			   			}}
			   			//}

			// cancer division
			   else if ((vacantSitesCancer(rI,rJ)>0) && (cell>0))
			   if (proliferation[cell]>=random.nextFloat()) {// If tossing the coin we are to proliferate...
				if ((cell==1) || ((cell==2) && (Age[rI][rJ]<maxProDivisions))) {// AND the cell is stem or TAC ...
				    //if (Oxygen[rI][rJ]>prolifThreshold) { // AND the oxygen concentration is enough for division..
					//consumption[rI][rJ]=consumptionDivision[Cells[rI][rJ]];
                    int[] daughter = findEmptySiteCancer (rI,rJ); // and there is space (for cancer)
			    	//    if ((daughter[0]==0) || (daughter[0]==size) || (daughter[1]==0)||(daughter[1]==size)) {simulationFinished=true;} // stop sim if a cell hits the edge
					births++; 
						if (cell==1) { // stem cell
							stemBirthsTotal[rI][rJ]++;
							stemBirthCounter[rI][rJ]++;
                            if (asymmetricRatio>random.nextFloat()) {
                            	Cells[daughter[0]][daughter[1]]=1;  // placing the stem daughter
                            	stemBirthCounter[daughter[0]][daughter[1]]=stemBirthCounter[rI][rJ]; //update stem birth counter
                            	carriedmutation[daughter[0]][daughter[1]]=carriedmutation[rI][rJ]; //inherit mutational status of parent
                            	carriedGenome[daughter[0]][daughter[1]]=carriedGenome[rI][rJ]; //inherit mutational status of parent
                            		if (mutfreq>random.nextFloat()) { // small chance of mutation
                            			mutationNum++; // advance mutation number
                            			System.out.println (+carriedmutation[rI][rJ]+", "+mutationNum+", "+timestep); // print (parent,child) pair
                            			//tree.put(carriedmutation[rI][rJ], mutationNum);
                            			//timeTree.put(mutationNum, timestep); // hash table stuff
                            			if (0.5>random.nextFloat()) {
                            				carriedmutation[daughter[0]][daughter[1]] = mutationNum;
                            				genomeToMod = new StringBuilder (carriedGenome[daughter[0]][daughter[1]]);
                            				genomeToMod.setCharAt(mutationNum-1,'1'); //  daughter carries new carriedGenome
                            				carriedGenome[daughter[0]][daughter[1]] = genomeToMod;
                            				//System.out.println (carriedGenome[rI][rJ]);
                            				//System.out.println (carriedGenome[daughter[0]][daughter[1]]);
                            				} // 50:50 mutate new position daughter
                            			else {
                            				carriedmutation[rI][rJ] = mutationNum;// else mutate original position daughter
                            				genomeToMod = new StringBuilder (carriedGenome[rI][rJ]);
                            				//genomeToMod = carriedGenome[rI][rJ];
                            				carriedGenome[rI][rJ]=genomeToMod;
                            				genomeToMod.setCharAt(mutationNum-1,'1');//  original carries new carriedGenome
                            				//System.out.println (carriedGenome[rI][rJ]);
                            				//System.out.println (carriedGenome[daughter[0]][daughter[1]]);
                            				}}}
                            else {Cells[daughter[0]][daughter[1]]=2; // asymmetric division, daughter is TAC
                            	stemBirthCounter[daughter[0]][daughter[1]]=0; // reset stem counter
                            	carriedmutation[daughter[0]][daughter[1]]=carriedmutation[rI][rJ]; // TAC carries parental mutation flag
                            	carriedGenome[daughter[0]][daughter[1]]=carriedGenome[rI][rJ];} // TAC carries parental genome
                                //} // redundant from above
				
								//else { // Only if there's hypoxia induced change of symmetric division ratio
                                //    float newASR=asymmetricRatio+0*(hypoxia-Oxygen[rI][rJ]); 
				    			// currently OFF by way of 0 the ^^ multiplier. 
                                //    if (newASR>random.nextFloat()) {Cells[daughter[0]][daughter[1]]=1;stemBirthCounter[daughter[0]][daughter[1]]=stemBirthCounter[rI][rJ];}
                                //    else {Cells[daughter[0]][daughter[1]]=2;stemBirthCounter[daughter[0]][daughter[1]]=0;} // Otherwise differentiate
                                //}
										} 
							else if (cell==2) {  // non-stem division
								//TACBirthCounter[rI][rJ]++;
								if (Age[rI][rJ]<maxProDivisions-1) {
									Cells[daughter[0]][daughter[1]]=2;
									Age[rI][rJ]++;
									Age[daughter[0]][daughter[1]]=Age[rI][rJ];
									carriedmutation[daughter[0]][daughter[1]]=carriedmutation[rI][rJ]; // TAC carries parental mutation flag
									carriedGenome[daughter[0]][daughter[1]]=carriedGenome[rI][rJ]; // TAC carries parental genome
								} else {
									Cells[daughter[0]][daughter[1]]=3;
									Cells[rI][rJ]=3;
									Age[rI][rJ]=0;
									Age[daughter[0]][daughter[1]]=Age[rI][rJ];
									carriedmutation[daughter[0]][daughter[1]]=carriedmutation[rI][rJ]; // TAC carries parental mutation flag
									carriedGenome[daughter[0]][daughter[1]]=carriedGenome[rI][rJ]; // TAC carries parental genome
								}
							}
						}
				} else if (pMotility>random.nextFloat()) { // Migration = not in use
						int[] daughter = findEmptySite (rI,rJ);
						Cells[daughter[0]][daughter[1]]=cell;
						Cells[rI][rJ]=0;
						Age[daughter[0]][daughter[1]]=Age[rI][rJ];
						Age[rI][rJ]=0;
						System.err.println ("moving "+rI+", "+rJ);
				}
			// Aging for mature cells
			if (cell==3) Age[rI][rJ]++;
	  	}
 	return true;
 }
    // this method finds vacant sites for HEALTHY cells
final int vacantSites (int x, int y) 

{
	// We assume that necrotic material counts as neighbour? NOTE: MOORE NEIGHBORHOOD
	int total=0;
	int[] p = new int [2];
		
	p=convertCoordinates (x+1,y-1);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x+1,y);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x+1,y+1);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x,y-1);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x,y+1);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x-1,y-1);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x-1,y);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x-1,y+1);
	if (Cells[p[0]][p[1]]==4) total++;
	return total;
}	

int[] findEmptySite (int x, int y)
{
	LinkedList vacantSites = new LinkedList();
	int[] tp1 = new int[2];
	int[] tp2 = new int[2];
	int[] tp3 = new int[2];
	int[] tp4 = new int[2];
	int[] tp5 = new int[2];
	int[] tp6 = new int[2];
	int[] tp7 = new int[2];
	int[] tp8 = new int[2];
	
	tp1=convertCoordinates (x+1,y-1);
	if (Cells[tp1[0]][tp1[1]]==4) vacantSites.add(tp1);
	tp2=convertCoordinates (x+1,y);
	if (Cells[tp2[0]][tp2[1]]==4) vacantSites.add(tp2);
	tp3=convertCoordinates (x+1,y+1);
	if (Cells[tp3[0]][tp3[1]]==4) vacantSites.add(tp3);
	tp4=convertCoordinates (x,y-1);
	if (Cells[tp4[0]][tp4[1]]==4) vacantSites.add(tp4);
	tp5=convertCoordinates (x,y+1);
	if (Cells[tp5[0]][tp5[1]]==4) vacantSites.add(tp5);
	tp6=convertCoordinates (x-1,y-1);
	if (Cells[tp6[0]][tp6[1]]==4) vacantSites.add(tp6);
	tp7=convertCoordinates (x-1,y);
	if (Cells[tp7[0]][tp7[1]]==4) vacantSites.add(tp7);
	tp8=convertCoordinates (x-1,y+1);
	if (Cells[tp8[0]][tp8[1]]==4) vacantSites.add(tp8);
	
	// Now let's see where.
	if (vacantSites.size() > 0) { // Now choose a vacant one, otherwise return the original location
		// pick a vacant site and return it
		int vacantElemIndex = random.nextInt(vacantSites.size());
		int[] p = (int[])vacantSites.get(vacantElemIndex);
		return (int[])p;	
	} else {
		int[] p = new int[2];
		p[0] = x; p[1] = y; // Just return the original
		System.out.println ("wrong!:"+vacantSites (x,y)+" - "+vacantSites.size());
		return p;
	}

}

    // a method for finding vacant sites for CANCER CELLS
final int vacantSitesCancer (int x, int y)
{
	// We assume that necrotic material and healthy cells DO NOT count as filling neighbouring sites. NOTE: MOORE NEIGHBORHOOD
	int total=0;
	int[] p = new int [2];
		
	p=convertCoordinatesNoFlux (x+1,y-1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinatesNoFlux (x+1,y);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinatesNoFlux (x+1,y+1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinatesNoFlux (x,y-1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinatesNoFlux (x,y+1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinatesNoFlux (x-1,y-1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinatesNoFlux (x-1,y);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinatesNoFlux (x-1,y+1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	return total;
	}

int[] findEmptySiteCancer (int x, int y)
{
	LinkedList vacantSitesCancer = new LinkedList();
	int[] tp1 = new int[2];
	int[] tp2 = new int[2];
	int[] tp3 = new int[2];
	int[] tp4 = new int[2];
	int[] tp5 = new int[2];
	int[] tp6 = new int[2];
	int[] tp7 = new int[2];
	int[] tp8 = new int[2];
	
	tp1=convertCoordinatesNoFlux (x+1,y-1);
	if ((Cells[tp1[0]][tp1[1]]==0) || (Cells[tp1[0]][tp1[1]]==4)) vacantSitesCancer.add(tp1);
	tp2=convertCoordinatesNoFlux (x+1,y);
	if ((Cells[tp2[0]][tp2[1]]==0) || (Cells[tp2[0]][tp2[1]]==4)) vacantSitesCancer.add(tp2);
	tp3=convertCoordinatesNoFlux (x+1,y+1);
	if ((Cells[tp3[0]][tp3[1]]==0) || (Cells[tp3[0]][tp3[1]]==4)) vacantSitesCancer.add(tp3);
	tp4=convertCoordinatesNoFlux (x,y-1);
	if ((Cells[tp4[0]][tp4[1]]==0) || (Cells[tp4[0]][tp4[1]]==4)) vacantSitesCancer.add(tp4);
	tp5=convertCoordinatesNoFlux (x,y+1);
	if ((Cells[tp5[0]][tp5[1]]==0) || (Cells[tp5[0]][tp5[1]]==4)) vacantSitesCancer.add(tp5);
	tp6=convertCoordinatesNoFlux (x-1,y-1);
	if ((Cells[tp6[0]][tp6[1]]==0) || (Cells[tp6[0]][tp6[1]]==4)) vacantSitesCancer.add(tp6);
	tp7=convertCoordinatesNoFlux (x-1,y);
	if ((Cells[tp7[0]][tp7[1]]==0) || (Cells[tp7[0]][tp7[1]]==4)) vacantSitesCancer.add(tp7);
	tp8=convertCoordinatesNoFlux (x-1,y+1);
	if ((Cells[tp8[0]][tp8[1]]==0) || (Cells[tp8[0]][tp8[1]]==4)) vacantSitesCancer.add(tp8);
	
	// Now let's see where.
	if (vacantSitesCancer.size() > 0) { // Now choose a vacant one, otherwise return the original location
		// pick a vacant site and return it
		int vacantElemIndexCancer = random.nextInt(vacantSitesCancer.size());
		int[] p = (int[])vacantSitesCancer.get(vacantElemIndexCancer);
		return (int[])p;	
		} else {
		int[] p = new int[2];
		p[0] = x; p[1] = y; // Just return the original
		System.out.println ("wrong!:"+vacantSitesCancer (x,y)+" - "+vacantSitesCancer.size());
		return p;
				}

}


    //////////////******************
/*
    public void iterateOxygen()
    {
        //
        float[][] newOxygen = new float[size][size];
        for (int rI = 0; rI < size; rI++)
            for (int rJ = 0; rJ < size; rJ++) {
                // Determine the actual coordinates for top (-1,0), left(0,-1), right(0,1), below(1,0)
                // using No Flux conditions
                int[] top = convertCoordinatesNoFlux(rI - 1, rJ);
                int[] left = convertCoordinatesNoFlux(rI, rJ - 1);
                int[] right = convertCoordinatesNoFlux(rI, rJ + 1);
                int[] below = convertCoordinatesNoFlux(rI + 1, rJ);

               // Diffusions
		if (Vasculature[rI][rJ] == 1) newOxygen[rI][rJ]=1;
		else { // do not update the vessels - they always stay 1
                newOxygen[rI][rJ]
                    = Oxygen[rI][rJ] + (kDe *
                                    (Oxygen[top[0]][top[1]]
                 + Oxygen[left[0]][left[1]]
                 + Oxygen[right[0]][right[1]]
                 + Oxygen[below[0]][below[1]]
                 - 4.0f * Oxygen[rI][rJ]));

                // Consumption
                if (Cells[rI][rJ]<4) { // if a cell is present
		    newOxygen[rI][rJ] = newOxygen[rI][rJ] - (consumption[rI][rJ])*((newOxygen[rI][rJ])/(newOxygen[rI][rJ] + Km)); ///Michaelis-Menten kinetics;  
                					}					
			}	 
			
            									}
        Oxygen = newOxygen;
    }
*/

//	public float [][] getOxygen() { return Oxygen; }
	public int [][] getCells () { return Cells; }
	//public int [][] getVasculature() {return Vasculature;}

    ///////  this is to run CA without Vis ***********  end from non sensitivty analysis stuff
     public static void main(String args[]) {
	 int maxTS=1000;
	      CA3 ca;
			System.err.println ("# CA version:"+CA3.version);
			float SCSymmetricDivBalance=0.2f;
				int maxDivisions=3;
				float densityV=0.04f;
					if (args.length==4) {
				    	SCSymmetricDivBalance = Float.parseFloat (args[0]);
						maxDivisions = Integer.parseInt (args[1]);
							maxTS=Integer.parseInt (args[2]);
							densityV=Float.parseFloat(args[3]);
								System.err.println ("Balance: "+SCSymmetricDivBalance+" maxDiv: "+maxDivisions+" maxTS: "+ maxTS);
									} else {
				    	System.err.println ("Arguments needed: s/a maxDivisions timesteps, densityV");
					System.exit(-1);
									}
				 ca = new CA3(SCSymmetricDivBalance, maxDivisions, densityV);
			 for (int ts=0;ts<maxTS;ts++) ca.nextTimeStep();
		    }
    ///////  this is to run CA without Vis ***********



/*
    ///////  this is to run CA without Vis ***********
    public static void main(String args[]) {
	 int maxTS=101;
	      CA3 ca;
			System.err.println ("# CA version:"+CA3.version);
			float asymmetricRatio=1.0f; // SCSymmetricDivBalance
			int maxProDivisions=10; // maxDivisions
			float densityVasculature=0.04f; // densityV
			Float DiffusionCoeff=0.1f;
			Float ConsumptionRate=1.0f;
			Float MichMen=0.1f;
			Float Hypox=0.1f;

				if (args.length==4) {
			    	DiffusionCoeff = Float.parseFloat (args[0]);
					ConsumptionRate = Float.parseFloat (args[1]);
					MichMen = Float.parseFloat(args[2]);
					Hypox = Float.parseFloat(args[3]);
					System.err.println ("Diffusion: "+DiffusionCoeff+" Consumption: "+ConsumptionRate+" MMconst: "+MichMen+" Hypoxia: "+Hypox);
						} else {
				    	System.err.println ("Arguments needed: diffusion consumption MichMen, hypoxia");
						System.exit(-1);
							}
				 ca = new CA3(DiffusionCoeff, ConsumptionRate, MichMen, Hypox);
			 for (int ts=0;ts<maxTS;ts++) ca.nextTimeStep();
		    }
    ///////  this is to run CA without Vis ***********
*/	

};
