/******************************************************************************

A Genetic Algorithm applies concepts (namely survival of the fittest and 
natural selection) from evolutionary theory to optimization problems. The 
Genetic Algorithm starts with a population of coded solutions (ChromosomePool) 
and evolves this population using the processes of Selection, Crossover and 
Mutation such that each successive genration of solutions is an improvement 
(on average) over previous generations.

******************************************************************************/
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "GeneticAlgorithm2.h"
#include "LatinHypercube.h"
#include "ParameterGroup.h"
#include <ParameterABC.h>
//#include "StatsClass.h"

#include "Exception.h"
#include "Utility.h"

#include <iostream>


/******************************************************************************
WarmStart()

Read the best solution from a previous run.
******************************************************************************/
/*void GeneticAlgorithm::WarmStart(void) {
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   double * pbest = new double[np+1];
   int newcount = SimpleWarmStart(np, pbest);
   m_pPopulation->SetChromosome(0, pbest);
   ((Model *)m_pModel)->SetCounter(newcount);
   delete [] pbest;
} */

/*
*****************************************************************************
CTOR

Initializes the population. First, all parameters are assigned default values
and then the user input file is checked for overriding values.
*****************************************************************************
*/
GeneticAlgorithm::GeneticAlgorithm() {

    FILE* pFile;
    double* pVals;
    char* pTok;
    LatinHypercube* pLHS = NULL;
    char* line;
    char tmp[DEF_STR_SZ];
    char tmp2[DEF_STR_SZ];
    int i, j, k, lvl, idx;
    IroncladString name = GetInFileName();

    // Read genetic algorithm information from the OSTRICH input file
    pFile = fopen(name, "r");
    if (pFile != NULL) {
        // Identify the Genetic Algorithm block within the file
        if (CheckToken(pFile, "BeginGeneticAlg", name) == true) {
            // Identify the end of the algorithm parameter block
            FindToken(pFile, "EndGeneticAlg", name);
            rewind(pFile);

            // Import data starting from the beginning of the block
            FindToken(pFile, "BeginGeneticAlg", name);
            line = GetNxtDataLine(pFile, name);

            while (strstr(line, "EndGeneticAlg") == NULL) {
                
                if (strstr(line, "PopulationSize") != NULL) {
                    // Set the the population size
                    sscanf(line, "%s %d", tmp, &m_NumPopulation);
                }
                else if (strstr(line, "MutationRate") != NULL) {
                    // Set the muation rate
                    sscanf(line, "%s %lf", tmp, &m_MutationRate);
                }
                else if (strstr(line, "Survivors") != NULL) {
                    // Set the number of survivors
                    sscanf(line, "%s %d", tmp, &m_NumSurvivors);
                }
                else if (strstr(line, "CrossoverRate") != NULL) {
                    // Set the number of survivors
                    sscanf(line, "%s %lf", tmp, &m_CrossoverRate);
                }
                else if (strstr(line, "NumGenerations") != NULL) {
                    // Set the number of generations
                    sscanf(line, "%s %d", tmp, &m_NumGenerationsMaximum);
                }
                else if (strstr(line, "InitPopulationMethod") != NULL) {
                    // Set the initializaiton method
                    sscanf(line, "%s %s", tmp, tmp2);
                    MyStrLwr(tmp2);
                    if (strcmp(tmp2, "random") == 0) { m_InitType = RANDOM_INIT; }
                    else if (strcmp(tmp2, "quadtree") == 0) { m_InitType = QUAD_TREE_INIT; }
                    else if (strcmp(tmp2, "lhs") == 0) { m_InitType = LHS_INIT; }
                }
                else if (strstr(line, "ConvergenceVal") != NULL) {
                    // Set the convergence value
                    sscanf(line, "%s %lf", tmp, &m_StopVal);
                }

                // Get the next line
                line = GetNxtDataLine(pFile, name);
            }
        }
        else {
            // The genetic algorithm has been specified but there is no input block. Log an error.
            LogError(ERR_FILE_IO, "Using default algorithm setup.");
        }

        /* initialize some or all pop. members to specied values */
        //rewind(pFile);
        //if (CheckToken(pFile, "BeginInitParams", name) == true) {
        //    FindToken(pFile, "EndInitParams", name);
        //    rewind(pFile);

        //    //allocate space for the parameter list
        //    num = m_pComm->GetParamGroupPtr()->GetNumParams();

        //    //count the number of entries
        //    FindToken(pFile, "BeginInitParams", name);
        //    line = GetNxtDataLine(pFile, name);
        //    m_NumInit = 0;
        //    while (strstr(line, "EndInitParams") == NULL) {
        //        m_NumInit++;
        //        line = GetNxtDataLine(pFile, name);
        //    }/* end while() */

        //    //allocate space for entries
        //    if (m_NumInit > 0) {
        //        NEW_PRINT("double *", m_NumInit);
        //        m_pInit = new double* [m_NumInit];
        //        MEM_CHECK(m_pInit);
        //        for (i = 0; i < m_NumInit; i++) {
        //            NEW_PRINT("double", num);
        //            m_pInit[i] = new double[num];
        //            MEM_CHECK(m_pInit[i]);
        //        }
        //    }/* end if() */

        //    //read in entries
        //    rewind(pFile);
        //    FindToken(pFile, "BeginInitParams", name);
        //    line = GetNxtDataLine(pFile, name);
        //    i = 0;
        //    while (strstr(line, "EndInitParams") == NULL) {
        //        pTok = line;
        //        //extract values, one-by-one, making any necessary conversions
        //        for (k = 0; k < num; k++) {
        //            j = ExtractString(pTok, tmp);
        //            j = ValidateExtraction(j, k, num, "ChromosomePool::Initialize()");
        //            pTok += j;
        //            m_pInit[i][k] = m_pComm->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
        //        }/* end for() */
        //        i++;
        //        line = GetNxtDataLine(pFile, name);
        //    }/* end while() */
        //}/* end if() */

        fclose(pFile);
    }/* end if() */


    /* 
    Perform setup checks of the configuration file
    */ 
    // Check that the population size is greater than zero
    if (m_NumPopulation <= 0) {
        LogError(ERR_FILE_IO, "Invalid population size");
        ExitProgram(1);
    }
    
    // Check that the mutation rate is valid
    if ((m_MutationRate < 0.00) || (m_MutationRate > 1.00)) {
        LogError(ERR_FILE_IO, "Invalid mutation rate");
        ExitProgram(1);
    }

    // Check that the crossover rate is valid
    if ((m_CrossoverRate < 0.00) || (m_CrossoverRate > 1.00)) {
        LogError(ERR_FILE_IO, "Invalid crossover rate");
        ExitProgram(1);
    }

    // Check that the number of generations is greater than zero
    if (m_NumGenerationsMaximum <= 0) {
        LogError(ERR_FILE_IO, "Invalid number of generations");
        ExitProgram(1);
    }


    /*
    Create the genes for the algorithm
    */
    // Create the genes
    m_pGenes = new RealEncodedGene[m_pParamGroup->GetNumParams()];
    MEM_CHECK(m_pGenes);

    // Loop and fill the data
    for (int entryParameter = 0; entryParameter < m_pParamGroup->GetNumParams(); entryParameter++) {
        
        // Get the parameter from the list
        ParameterABC* test = m_pParamGroup->GetParamPtr(entryParameter);

        // Set the parameter values into the Gene 
        m_pGenes[entryParameter].SetValue(test->GetInitialValueTransformed());
        m_pGenes[entryParameter].SetUpr(test->GetUpperBoundTransformed());
        m_pGenes[entryParameter].SetLwr(test->GetLowerBoundTransformed());
        m_pGenes[entryParameter].SetMutationRate(m_MutationRate);
        m_pGenes[entryParameter].SetCrossoverRate(m_CrossoverRate);

    }

   
   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()

Free up memory used by the GA and it's member variables.
******************************************************************************/
void GeneticAlgorithm::Destroy(void) {
                                                 
                       
    delete m_pGenes;

                


   IncDtorCount();
}

/*
Initialize()

*/
double** GeneticAlgorithm::CreateInitialSample(int sampleSize) {

    // Create the holder for the samples
    int numberOfParameters = m_pParamGroup->GetNumParams();
    double** samples;
    
    // Generate the set of initial parameters
    if (m_InitType == LHS_INIT) {
        // Use Latin Hypercube to initalize the space
        // Create a new LHS object
        LatinHypercube *pLHS = new LatinHypercube(sampleSize, numberOfParameters, numberOfParameters);
        MEM_CHECK(pLHS);

        // Extract the bounds needed to create the LHS matrix
        double* lowerBounds = new double[numberOfParameters];
        double* upperBounds = new double[numberOfParameters];
        
        for (int i = 0; i < numberOfParameters; i++) {
            lowerBounds[i] = m_pGenes[i].GetLwr();
            upperBounds[i] = m_pGenes[i].GetUpr();
        }

        // Fill the array from the LHS object
        pLHS->CreateUniformSample(lowerBounds, upperBounds);
        samples = pLHS->GetSampleMatrix();
        
        // Cleanup the memory
        //delete pLHS;
        delete[] lowerBounds;
        delete[] upperBounds;

    }
    else if (m_InitType == RANDOM_INIT){
        // Use random values to initialize the space
        // Create the sample holder
        samples = new double* [sampleSize];

        // Loop on each chromosome/gene and draw a random value for each.
        for (int entryChromosome = 0; entryChromosome < sampleSize; entryChromosome++) {
            // Create the holder for this chromosome of samples
            samples[entryChromosome] = new double[numberOfParameters];

            for (int entryParameter = 0; entryParameter < numberOfParameters; entryParameter++) {
                samples[entryChromosome][entryParameter] = m_pGenes[entryParameter].GetRandomValue();
            }
        }
    }
    //else if (m_InitType == QUAD_TREE_INIT) {
    //    // Use quad trees to initialize the space
    //    // Create the sample holder
    //    samples = new double* [sampleSize];

    //    // Create the new QuadTree object
    //    QuadTree *m_pTrees = new QuadTree[numberOfParameters];

    //    // Set the tree values
    //    int lvl = 0;
    //    int idx = 0;

    //    // Loop over the tree and generate values
    //    for (int entryChromosome = 0; entryChromosome < sampleSize; entryChromosome++) {
    //        // Create the holder for this chromosome of samples
    //        samples[entryChromosome] = new double[numberOfParameters];

    //        // Fill the data into the trees
    //        for (int j = 0; j < numberOfParameters; j++) {
    //            double lwr = m_pGenes[j].GetLwr();
    //            double upr = m_pGenes[j].GetUpr();
    //            m_pTrees[j].Init(lwr, upr);
    //        }

    //        // Call the head of the tree
    //        double* pVals = GetTreeCombo(lvl, idx, m_pTrees, numberOfParameters);

    //        // Contine to loop on the tree structure
    //        if (pVals == NULL) {
    //            for (int j = 0; j < numberOfParameters; j++) { m_pTrees[j].Expand(); }
    //            lvl++;
    //            idx = 0;
    //            pVals = GetTreeCombo(lvl, idx, m_pTrees, numberOfParameters);
    //        }

    //        idx++;

    //        samples[entryChromosome] = pVals;
    //        delete[] pVals;
    //    }
    //}
    else {
        // Log an error
         // TODO: Add this
        samples = NULL;
    }
 
    // Return to the calling function
    return samples;
    
}




/*
*****************************************************************************
TourneySelection()

Determines the mating pool by randomly selecting two chromosomes and comparing
their fitness values. The chromosome with the better fitness gains the right
to pass its genes into the next generation. The configuration variable
m_NumSurvivors is used to guarantee that the top chromosomes survive unchanged
into the next generation.

The input argument specifies the number of combatants in the tournament. For
'standard' GA this is set equal to 2. For computation constrained, the number
of combatants increases as number of generations increases.
*****************************************************************************
*/
void GeneticAlgorithm::TourneySelection(int nCombatants, double* objectives, int numberOfObjectives, double** samples, double** scratch) {
    int r1, r2;
    double* pMax;   //chromosome with max. fitness
    double fit1;
    double fit2;
    double maxFit;
    double lastMax;


    /*-------------------------------------------
    Reserve the top m_NumSurvivors chromosomes.
    --------------------------------------------*/
    lastMax = NEARLY_HUGE;
    pMax = NULL;
    for (int i = 0; i < m_NumSurvivors; i++) {
        maxFit = -NEARLY_HUGE;

        for (int j = 0; j < m_NumPopulation; j++) {
            fit1 = objectives[j];

            if ((fit1 > maxFit) && (fit1 <= lastMax) && (pMax != scratch[i])) {
                pMax = samples[j];
                maxFit = fit1;
            }
        }

        //propagate nth max. to next generation
        lastMax = maxFit;
        scratch[i] = pMax;
    }

    /*-------------------------------------------
    Use n-member tourney to select the remaining chromosomes.
    --------------------------------------------*/
    for (int i = m_NumSurvivors; i < m_NumPopulation; i++) {
        // pick random chromosomes 
        r1 = MyRand() % m_NumPopulation;
        fit1 = objectives[r1];

        for (int j = 0; j < (nCombatants - 1); j++) {
            r2 = MyRand() % m_NumPopulation;

            //evaluate their fitness
            fit2 = objectives[r1];

            //the better one gets to go to the nextGeneration
            if (fit2 > fit1) {
                fit1 = fit2;
            }
        }

        scratch[i] = samples[r1];
    } 
} 

/*
*****************************************************************************
Crossover()

Crosses over each chromsome of the population with the next one in the
population, except those that are in the top <m_NumSurvivors>.
*****************************************************************************
*/
void GeneticAlgorithm::Crossover(double* objectives, int numberOfObjectives, double** samples, double** scratch) {

    // Add the best back to the worker to maintain it in the population
    scratch[0] = m_BestAlternative;

    // Find the best from the previous generation to maintain in the current generation
    if (m_NumSurvivors > 1) {
        // Create the arrays
        double *objectiveValues = new double[m_NumSurvivors - 1];
        int *objectiveIndices = new int[m_NumSurvivors - 1];

        // Fill them with values
        for (int entryIndex = 0; entryIndex < m_NumSurvivors - 1; entryIndex++) {
            objectiveValues[entryIndex] = INFINITY;
            objectiveIndices[entryIndex] = -1;
        }

        for (int entryObjective = 0; entryObjective < numberOfObjectives; entryObjective++) {
            double maxDifference = 0;
            int maxDifferenceIndex = -1;

            for (int entryIndex = 0; entryIndex < m_NumSurvivors - 1; entryIndex++) {
                double difference = objectiveValues[entryIndex] - objectives[entryObjective];

                if (difference > maxDifference) {
                    maxDifference = difference;
                    maxDifferenceIndex = entryIndex;
                }
            }

            objectiveIndices[maxDifferenceIndex] = maxDifferenceIndex;
            objectiveValues[maxDifferenceIndex] = objectives[entryObjective];

        }

        // Set the remaining survivors into the scratch array
        for (int entryIndex = 1; entryIndex < m_NumSurvivors - 1; entryIndex++) {
            scratch[entryIndex] = samples[objectiveIndices[entryIndex - 1]];
        }
    }

    // Adjust for keepting the best in the pool if there are no survivors specified
    int survivorsAdjusted = 1;
    if (m_NumSurvivors > survivorsAdjusted) {
        survivorsAdjusted = m_NumSurvivors;
    }
    
    //crossover everyone with their neighbor
    for (int i = survivorsAdjusted; i < numberOfObjectives - 1; i++) {
        // Get the parent chromosomes
        double* pMom = samples[i];
        double* pPop = samples[i + 1];

        // Get the parent objectives
        double pMomObjective = objectives[i];
        double pPopObjective = objectives[i + 1];

        // Create the array to hold the new chromosome
        double* tempChromosome = new double[m_pParamGroup->GetNumParams()];

        // Crossover each entry in the array
        for (int entryParam = 0; entryParam < m_pParamGroup->GetNumParams(); entryParam++) {
            tempChromosome[entryParam] = m_pGenes[entryParam].Crossover(&m_pGenes[entryParam + 1], samples[i][entryParam], samples[i + 1][entryParam], 
                                                                         pMomObjective, pPopObjective, m_pParamGroup->GetNumParams());

        }

        // Set the chromosome into the scratch array
        scratch[i] = tempChromosome;
    }
}


/*
*****************************************************************************
Mutate()

Mutates individual chromsomes of the population according to a
pre-established mutation rate.
*****************************************************************************
*/
void GeneticAlgorithm::Mutate(double** scratch, int numberOfSamples) {


    for (int entryAlternative = m_NumSurvivors; entryAlternative < m_NumPopulation; entryAlternative++) {
        for (int entryParam = 0; entryParam < m_pParamGroup->GetNumParams(); entryParam++) {
            scratch[entryAlternative][entryParam] = m_pGenes[entryParam].Mutate(scratch[entryAlternative][entryParam]);
        }
    }
}


/*
*****************************************************************************
CreateSample()

Creates the next generation of the chromosome population using tourney selection,
crossover, and mutation.
*****************************************************************************
*/
double** GeneticAlgorithm::CreateSample(double* objectives, int numberOfObjectives, double** samples) {

    // Create a scratch array to hold the updated dated
    double** scratch; 
    scratch = new double* [numberOfObjectives];
    for (int entrySample = 0; entrySample < numberOfObjectives; entrySample++) {
        scratch[entrySample] = new double[numberOfObjectives];
    }
        
    // Go through the genetic sampling process
    TourneySelection(2, objectives, numberOfObjectives, samples, scratch);
    Crossover(objectives, numberOfObjectives,  samples, scratch);
    Mutate( scratch, numberOfObjectives);

    // Swap scratch with the samples
    return scratch;
}


/*
*****************************************************************************
FreezeGenes()

For each chromosome, randomly freeze the given number of genes so that they
are at the current global optimal.
*****************************************************************************
*/
/*void GeneticAlgorithm::FreezeGenes(int numFreeze) {
    int np = m_pComm->GetParamGroupPtr()->GetNumParams();

    for (int i = m_NumSurvivors; i < m_PoolSize; i++) {
        SampleWithReplacement(-2, np);
        for (int j = 0; j < numFreeze; j++) {
            int k = SampleWithReplacement(1, np);
            m_pPool[i]->GetGenePtr(k)->SetValue(GetBestFit()->GetGenePtr(k)->GetValue());
        }
    }
}*/

/*
*****************************************************************************
CalcAvgFitness()

Calculates and returns the average fitness of the population. This parameter
is used in the termination criteria of the genetic algorithm.
*****************************************************************************
*/
double GeneticAlgorithm::CalcMeanFitness(double* objectives, int numberOfObjectives) {

    // Define the holder for the sum
    double sum = 0.00;

    // Loop over each entry and sum the value
    for (int i = 0; i < numberOfObjectives; i++) { 
        sum += objectives[i];
    }

    // Divide the total by the number of objectives to get the value
    double avg = sum / (double)numberOfObjectives;

    return avg;
} 

/*
*****************************************************************************
CalcMedianFitness()

Calculates and returns the median fitness of the population. This parameter
is used in the termination criteria of the genetic algorithm.
*****************************************************************************
*/
double GeneticAlgorithm::CalcMedianFitness(double* objectives, int numberOfObjectives) {
    
    double med = CalcMedian(objectives, numberOfObjectives);

    return med;
} 

/*
*****************************************************************************
GetBestFit()

Retrieves the chromosome value and index that has the best fitness value.
*****************************************************************************
*/
void GeneticAlgorithm::GetBestObjective(double* objectives, int numberOfObjectives) {
    
    double iterationMinimum = objectives[0];
    int iterationIndex = 0;

    for (int i = 1; i < numberOfObjectives; i++) {
        if (objectives[i] < iterationMinimum) {
            iterationMinimum = objectives[i];
            iterationIndex = i;
        }
    }

    m_BestObjectiveIteration = iterationMinimum;
    m_BestObjectiveIndexIteration = iterationIndex;

}


/******************************************************************************
EvalFitSuperMUSE()

Compute fitness of entire population using SuperMUSE. This routine interfaces
with the RepeatTasker SuperMUSE program, which assigns model evaluations to
SuperMUSE clients on a first-come-first-served basis.
******************************************************************************/
//void ChromosomePool::EvalFitSuperMUSE(void) {
//    double val;
//    bool bOk;
//    int pop_size;
//    int i;
//    ParameterGroup* pGroup;
//    SuperMUSE* pSMUSE = GetSuperMusePtr();
//
//    pop_size = m_PoolSize;
//
//    /* ----------------------------------------------------------------
//    Generate task file that describes the desired parallel evaluations.
//    This is somewhat analogous to the BcastPopulation() operation used
//    for MPI-parallel operations. Write the parameter values of each
//    population member as entries in the task file.
//
//    Entries are first accumlated into a temp file to prevent the
//    SuperMUSE RepeatTasker program from prematurely processing the task
//    file.
//    ---------------------------------------------------------------- */
//    for (i = 0; i < pop_size; i++)
//    {
//        //stuff the parameter group with values
//        pGroup = ConvertChromosome(m_pPool[i]);
//
//        //pass group to supermuse
//        pSMUSE->WriteTask(pGroup);
//    }/* end for() */
//
//    //Finish task file (this will cause RepeatTasker to begin processing the job)
//    pSMUSE->FinishTaskFile();
//
//    //wait for SuperMUSE to report back (via the success or error files)
//    bOk = pSMUSE->WaitForTasker();
//
//    if (bOk == false) //SuperMUSE failed
//    {
//        LogError(ERR_SMUSE, "Reverting to serial execution.");
//        DisableSuperMUSE();
//        EvalFitness();
//    }
//    else //SuperMUSE was successful
//    {
//        for (i = 0; i < pop_size; i++)
//        {
//            /* -----------------------------------------------
//            Stuff the parameter group with ith population
//            member. This ensures that each objective function
//            gets associated with the correct parameter values.
//            ------------------------------------------------ */
//            pGroup = ConvertChromosome(m_pPool[i]);
//
//            //stuff i-th result into chromosome pool
//            val = pSMUSE->GatherResult(i);
//            m_pPool[i]->SetFitness(-val);
//        }/* end for() */
//    }/* end else() */
//}/* end EvalFitSuperMUSE() */



/******************************************************************************
Initialize()

Initializes the population so that it is on a budget.
******************************************************************************/
//void ChromosomePool::Initialize(int* budget)
//{
//    FILE* pFile;
//    int popSize, num, np;
//    double rate, lwr, upr;
//    double* pVals;
//    char* pTok;
//    LatinHypercube* pLHS = NULL;
//    char* line;
//    char tmp[DEF_STR_SZ];
//    int i, j, k, lvl, idx;
//    IroncladString name = GetInFileName();
//
//    //assign default GA parameters
//    np = m_pComm->GetParamGroupPtr()->GetNumParams();
//    *budget = 1000;
//    m_NumGenerations = (int)(0.5 + 2 * np + sqrt((double)*budget));
//    m_NumInit = 0;
//    popSize = (int)(0.5 + (*budget / m_NumGenerations));
//    rate = 0.15; //initial mutation rate is high (15%)
//    m_NumSurvivors = (int)(MyMax(1.00, 0.5 + 0.05 * popSize));
//    m_InitType = LHS_INIT;
//    m_StopVal = -1.00;
//
//    pFile = fopen(name, "r");
//    if (pFile != NULL)
//    {
//        if (CheckToken(pFile, "BeginGeneticAlg", name) == true)
//        {
//            FindToken(pFile, "EndGeneticAlg", name);
//            rewind(pFile);
//
//            FindToken(pFile, "BeginGeneticAlg", name);
//            line = GetNxtDataLine(pFile, name);
//
//            while (strstr(line, "EndGeneticAlg") == NULL)
//            {
//                if (strstr(line, "Budget") != NULL)
//                {
//                    sscanf(line, "%s %d", tmp, budget);
//                    if (*budget <= 0) *budget = 1000;
//                }
//                line = GetNxtDataLine(pFile, name);
//            }/* end while() */
//        }/* end if() */
//        else
//        {
//            LogError(ERR_FILE_IO, "Using default algorithm setup.");
//        }/* end else() */
//
//        /* initialize some or all pop. members to specied values */
//        rewind(pFile);
//        if (CheckToken(pFile, "BeginInitParams", name) == true)
//        {
//            FindToken(pFile, "EndInitParams", name);
//            rewind(pFile);
//
//            //allocate space for the parameter list
//            num = m_pComm->GetParamGroupPtr()->GetNumParams();
//
//            //count the number of entries
//            FindToken(pFile, "BeginInitParams", name);
//            line = GetNxtDataLine(pFile, name);
//            m_NumInit = 0;
//            while (strstr(line, "EndInitParams") == NULL)
//            {
//                m_NumInit++;
//                line = GetNxtDataLine(pFile, name);
//            }/* end while() */
//
//            //allocate space for entries
//            if (m_NumInit > 0)
//            {
//                NEW_PRINT("double *", m_NumInit);
//                m_pInit = new double* [m_NumInit];
//                MEM_CHECK(m_pInit);
//                for (i = 0; i < m_NumInit; i++)
//                {
//                    NEW_PRINT("double", num);
//                    m_pInit[i] = new double[num];
//                    MEM_CHECK(m_pInit[i]);
//                }
//            }/* end if() */
//
//            //read in entries
//            rewind(pFile);
//            FindToken(pFile, "BeginInitParams", name);
//            line = GetNxtDataLine(pFile, name);
//            i = 0;
//            while (strstr(line, "EndInitParams") == NULL)
//            {
//                pTok = line;
//                //extract values, one-by-one, making any necessary conversions
//                for (k = 0; k < num; k++)
//                {
//                    j = ExtractString(pTok, tmp);
//                    j = ValidateExtraction(j, k, num, "ChromosomePool::Initialize()");
//                    pTok += j;
//                    m_pInit[i][k] = m_pComm->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
//                }/* end for() */
//                i++;
//                line = GetNxtDataLine(pFile, name);
//            }/* end while() */
//        }/* end if() */
//
//        fclose(pFile);
//    }/* end if() */
//
//    /* -------------------------------------------------------------------
//    Adjust population size and max. gens to reflect user-defined budget
//    ------------------------------------------------------------------- */
//    if (*budget > popSize * m_NumGenerations)
//    {
//        //inc. max. gens
//        m_NumGenerations = *budget / popSize;
//    }
//    else if (*budget < (popSize * m_NumGenerations))
//    {
//        if (*budget < (popSize * 3)) //revise pop. size
//        {
//            popSize = (int)MyMax(*budget / 3, 3);
//            m_NumGenerations = *budget / popSize;
//        }
//        else //revise max gens
//        {
//            m_NumGenerations = *budget / popSize;
//        }
//    }
//    if (popSize * m_NumGenerations < *budget) m_NumGenerations++;
//
//    m_Generation = 0;
//    m_Proto = m_pComm->CreateProto(rate);
//
//    //init. mutation counters
//    NEW_PRINT("int", m_Proto->GetNumGenes());
//    m_pMutCount = new int[m_Proto->GetNumGenes()];
//    MEM_CHECK(m_pMutCount);
//    for (i = 0; i < m_Proto->GetNumGenes(); i++) { m_pMutCount[i] = 0; }
//
//    m_PoolSize = popSize;
//    NEW_PRINT("Chromosome *", m_PoolSize);
//    m_pPool = new Chromosome * [m_PoolSize];
//    MEM_CHECK(m_pPool);
//
//    NEW_PRINT("Chromosome *", m_PoolSize);
//    m_pScratch = new Chromosome * [m_PoolSize];
//    MEM_CHECK(m_pScratch);
//
//    NEW_PRINT("double", m_PoolSize);
//    m_Fmedian = new double[m_PoolSize];
//    MEM_CHECK(m_Fmedian);
//
//    NEW_PRINT("double", m_Proto->GetNumGenes());
//    pVals = new double[m_Proto->GetNumGenes()];
//    MEM_CHECK(pVals);
//
//    NEW_PRINT("LatinHypercube", 1);
//    pLHS = new LatinHypercube(m_Proto->GetNumGenes(), m_PoolSize);
//    MEM_CHECK(pLHS);
//
//    for (j = 0; j < m_Proto->GetNumGenes(); j++)
//    {
//        lwr = m_Proto->GetGenePtr(j)->GetLwr();
//        upr = m_Proto->GetGenePtr(j)->GetUpr();
//        pLHS->InitRow(j, lwr, upr);
//    }/* end for() */
//
//    lvl = idx = 0;
//    for (i = 0; i < m_PoolSize; i++)
//    {
//        for (j = 0; j < m_Proto->GetNumGenes(); j++) { pVals[j] = pLHS->SampleRow(j); }
//        m_pPool[i] = m_Proto->CreateChromo(pVals);
//        m_pScratch[i] = m_Proto->CreateChromo(pVals);
//    }/* end for() */
//
//    //seed initial population
//    for (i = 0; i < m_NumInit; i++)
//    {
//        delete m_pPool[i];
//        m_pPool[i] = m_Proto->CreateChromo(m_pInit[i]);
//        delete m_pScratch[i];
//        m_pScratch[i] = m_Proto->CreateChromo(m_pInit[i]);
//    }
//    delete pLHS;
//    delete[] pVals;
//} /* end Initialize() */



/******************************************************************************
WriteMetrics()

Write out setup and metrics for the pool.
******************************************************************************/
void GeneticAlgorithm::WriteStartingMetrics(void) {
    // TODO: Move to the write utility class

    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Log the algorithm setup information
    fprintf(pFile, "Population Size         : %d\n", m_NumPopulation);
    fprintf(pFile, "Mutation Rate           : %f\n", m_MutationRate);
    fprintf(pFile, "Crossover Rate          : %f\n", m_CrossoverRate);
    fprintf(pFile, "Number of Elites        : %d\n", m_NumSurvivors);
    fprintf(pFile, "Initialization Method   : ");
    if (m_InitType == RANDOM_INIT) { fprintf(pFile, "Random\n"); }
    else if (m_InitType == QUAD_TREE_INIT) { fprintf(pFile, "Quad-Tree\n"); }
    else if (m_InitType == LHS_INIT) { fprintf(pFile, "Latin Hypercube Sampling\n"); }
    else { fprintf(pFile, "Unknown\n\n"); }

    // Write the parameter to the file
    fprintf(pFile, "Run   "); // TODO: align columns 
    fprintf(pFile, "%-12s  ", "Objective");
    m_pParamGroup->Write(pFile, WRITE_BNR);
    fprintf(pFile, "Convergence\n");    

    // Close the log file
    fclose(pFile);


}/* end WriteMetrics() */

/******************************************************************************
WriteEndingMetrics()

Write out setup and metrics for the algorithm.
******************************************************************************/
void GeneticAlgorithm::WriteEndingMetrics(void) {
    //TODO: Move htis to the write utility class
    
    // Open the log file
    char fileName[DEF_STR_SZ];
    FILE* pFile;
    sprintf(fileName, "OstOutput%d.txt", 0);
    pFile = fopen(fileName, "a");

    // Write the result information
    WriteOptimalToFileWithGroup2(pFile, m_pParamGroup, m_BestObjective);
    
    // Write the algorithm information
    fprintf(pFile, "\nAlgorithm Metrics\n");
    fprintf(pFile, "Algorithm               : GeneticAlgorithm\n");
    fprintf(pFile, "Desired Convergence Val : %E\n", m_StopVal);
    fprintf(pFile, "Actual Convergence Val  : %E\n", m_CurStop);
    fprintf(pFile, "Max Generations         : %d\n", m_NumGenerationsMaximum);
    fprintf(pFile, "Actual Generations      : %d\n", m_Generation);

    if (m_CurStop <= m_StopVal) {
        fprintf(pFile, "Algorithm successfully converged on a solution\n");
    }
    else {
        fprintf(pFile, "Algorithm failed to converge on a solution, more generations may be needed\n");
    }

    // Close the log file
    fclose(pFile);
}


/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using the GA.
******************************************************************************/
/*void GeneticAlgorithm::Calibrate(void) { 
   int id;
   char fileName[DEF_STR_SZ];
   FILE * pFile;

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   RegisterStatsPtr(m_pStats);

   Optimize();
   
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //compute statistics (variance and covariance)
   m_pStats->CalcStats();

   if(id == 0) {
      sprintf(fileName, "OstOutput%d.txt", id);

      //write statistics of best parameter set to output file
      pFile = fopen(fileName, "a");   
      m_pStats->WriteStats(pFile);
      fclose(pFile);

      //write statistics of best parameter set to output file
      m_pStats->WriteStats(stdout);
   }
}*/

/******************************************************************************
Optimize()

Minimize the objective function using the GA.
******************************************************************************/
void GeneticAlgorithm::Optimize(void) {

    // Write the information to the primary log file
    WriteSetup2(this, "GeneticAlgorithm");
    WriteStartingMetrics();

    // Initialize the workers
    ConfigureWorkers();

   // Initialize the chromosome for the first solve
   double** samples = CreateInitialSample(m_NumPopulation);

    // Enter the solution loop
   while (m_NumGenerationsMaximum > m_Generation && m_CurStop >= m_StopVal) {

       // Evaluate the fitness of the sample set
       double *objectives = new double[m_NumPopulation];
       ManageSingleObjectiveIterations(samples, m_pParamGroup->GetNumParams(), m_NumPopulation, objectives);

       // Compute the best and averages
       GetBestObjective(objectives, m_NumPopulation);
       double meanObjective = CalcMeanFitness(objectives, m_NumPopulation);
       double medianOBjective = CalcMedianFitness(objectives, m_NumPopulation);
       m_CurStop = fabs((medianOBjective - m_BestObjectiveIteration) / medianOBjective);

       // Swap the current with the best if doing better
       if (m_BestObjectiveIteration < m_BestObjective) {
           m_BestObjective = m_BestObjectiveIteration;
           m_BestAlternative = samples[m_BestObjectiveIndexIteration];
       }
       
       // Log the iteration
       // Set the values into the pointer groups
       for (int entryParam = 0; entryParam < m_pParamGroup->GetNumParams(); entryParam++) {
           ParameterABC* temp = m_pParamGroup->GetParamPtr(entryParam);
           temp->SetEstimatedValueTransformed(m_BestAlternative[entryParam]);
       }

       // Write iterationr result to the log file
       WriteRecord2(this, m_Generation, m_BestObjective, m_CurStop);

       // Increment the generation variable
       m_Generation++;

       // Generate another set of samples if solution will continue
       if (m_Generation < m_NumGenerationsMaximum && m_CurStop >= m_StopVal) {
           samples = CreateSample(objectives, m_NumPopulation, samples);
       }
   }


   // Terminate the secondary workers and let the clean up.
   TerminateWorkers();

   //write algorithm metrics
   WriteEndingMetrics();

   // Cleanup memory space
   delete samples;
} 

/******************************************************************************
GA_Program()

Calibrate the model using the GA.
******************************************************************************/
void GA_Program(int argC, StringType argV[]) {

   // Get the rank of the process
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   // Determine behavior based on the rank
   if (rank == 0) {
       // This is the primary worker. Setup the algorithm and call the processing routine
       // Create the genetice agorithm
       GeneticAlgorithm *GA = new GeneticAlgorithm();

       // Begin the primary worker operations. This transfers information to the secondary workers
       // and manages the solution process
       GA->Optimize();



   }
   else {
       // This is a secondary worker. Create the secondary worker and call the work routine
       // Create and setup the secondary worker
       ModelWorker worker = ModelWorker();

       // Start the work in the processor 
       worker.Work();

       // Tear-down is triggered by the primary worker
   }

} 


