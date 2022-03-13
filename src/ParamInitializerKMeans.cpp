/******************************************************************************
File      : KMeansParamInitializer.cpp
Author    : L. Shawn Matott
Copyright : 2022, L. Shawn Matott

Encapsulates the KMeans clustering parameter initialization strategy for use
with LID placement in the EPA Stormwater Management Model.

Version History
02-27-2022    lsm   created
******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include "ParamInitializerABC.h"
#include "ParameterABC.h"
#include "ParameterGroup.h"
#include "Exception.h"
#include "Utility.h"

/* ----------------------------------------------------------------------------
   GetClusterArea()

    Lookup the area of "cluster_id".
---------------------------------------------------------------------------- */
double KMeansParamInitializer::GetClusterArea(SubcatchmentCoordStruct * pData, int n, int cluster_id)
{
    int i;

    for(i = 0; i < n; i++)
    {
        if(pData[i].cluster_id == cluster_id)
        {
            return(pData[i].cluster_area_ac);
        }
    }
    return(0.00);
}/* end GetClusterArea() */

/* ----------------------------------------------------------------------------
   pick_random_cluster()

    Pick a random cluster with replacement.
---------------------------------------------------------------------------- */
int KMeansParamInitializer::PickRandomCluster(int n)
{
    FILE * pIn;
    int i;
    int r;
    char rstr[100];

    // for now read the n-th value from file
    pIn = fopen("random_clusters.out", "r");
    for(i = 0; i < n; i++)
    {
        fscanf(pIn, "%s\n", rstr);
        r = atoi(rstr);
    }
    fclose(pIn);
    return r;
}/* end PickRandomCluster() */

/* -------------------------------------------------------------------------------------------
   AppendClusterAreas()

    Associate cluster labels with total area of associated subcatchments.
-------------------------------------------------------------------------------------------- */
SubcatchmentCoordStruct * KMeansParamInitializer::
      AppendClusterAreas(SubcatchmentCoordStruct * data, int n_subcatchments, int n_clusters)
{
    int i;
    int j;
    double cluster_area;


    for(i = 0; i < n_clusters; i++)
    {
        // accumulate the area of the i-th cluster
        cluster_area = 0.00;
        for(j = 0; j < n_subcatchments; j++)
        {
            if(data[j].cluster_id == i)
            {
                cluster_area += data[j].area_ac;
            }
        }
        // assign the aggregated area to data points in the cluster
        for(j = 0; j < n_subcatchments; j++)
        {
            if(data[j].cluster_id == i)
            {
                data[j].cluster_area_ac = cluster_area;
            }
        }    
    }
    return data;
}/* end AppendClusterAreas() */

/* ----------------------------------------------------------------------------
   AppendLabels()

    Associate cluster labels with subcatchment data.
---------------------------------------------------------------------------- */
SubcatchmentCoordStruct * KMeansParamInitializer::
              AppendLabels(SubcatchmentCoordStruct * data, int * labels, int n)
{
    int i;

    for(i = 0; i < n; i++)
    {
        data[i].cluster_id = labels[i];
    }
    return data;
}/* end AppendLabels() */

/* ----------------------------------------------------------------------------
   CreateClusters()

    Analyze coordinate data and organize into clusters.
---------------------------------------------------------------------------- */
void KMeansParamInitializer::CreateClusters(SubcatchmentCoordStruct * data, KMeansClusterStruct * pKMeans)
{
    int i;
    int j;
    int r;
    int nData = pKMeans->num_data_points;
    int nClust = pKMeans->num_clusters;
    int count;
    double xc;
    double yc;
    FILE * pIn;

    // for now read assignments from file
    pIn = fopen("cluster_map_ipynb.out", "r");
    for(i = 0; i < nData; i++)
    {
        // r = (rand() % nClust) + 1;
        fscanf(pIn, "%d,%d", &j, &r);
        pKMeans->labels[i] = r;
    }
    fclose(pIn);
    // compute cluster centers
    for(i = 0; i < nClust; i++)
    {
        count = 0;
        xc = 0.00;
        yc = 0.00;

        // aggregate (x,y) values for data points in the ith cluster
        for(j = 0; j < nData; j++)
        {
            if(i == pKMeans->labels[j])
            {
                count++;
                xc += data[j].x;
                yc += data[j].y;
            }
        }
        // compute the center points usinmg average value
        if(count != 0)
        {
            pKMeans->centers[i].x = xc / (float) count;
            pKMeans->centers[i].y = yc / (float) count;
        }
        else
        {
            printf("WARNING - cluster %d has no members!", i);
            pKMeans->centers[i].x = 0.00;
            pKMeans->centers[i].y = 0.00;
        }
    }

    return;
}/* end CreateKmeansClusters() */

/* ----------------------------------------------------------------------------
   rtrim()

    Remove carriage return and line feed from end of a string.
---------------------------------------------------------------------------- */
void KMeansParamInitializer::rtrim(char * pStr)
{
    pStr[strcspn(pStr, "\r\n")] = 0x0;
}/* end rtrim() */

/* ----------------------------------------------------------------------------
   PrintSubcatchmentData()

    A function that displays subcatchment data
---------------------------------------------------------------------------- */
void KMeansParamInitializer::PrintSubcatchmentData(SubcatchmentCoordStruct * data, int n)
{
    fprintf(m_pLog, "%-30s\t%-14s\t%-14s\t%s\n","NAME","X","Y","AREA");
    if(n < 10)
    {
        for(int i = 0; i < n; i++)
        {
            fprintf(m_pLog, "%-30s\t%f\t%f\t%f\n", data[i].name, data[i].x, data[i].y, data[i].area_ac);
        }
    }
    else
    {
        // first 5
        for(int i = 0; i < 5; i++)
        {
            fprintf(m_pLog, "%-30s\t%f\t%f\t%f\n", data[i].name, data[i].x, data[i].y, data[i].area_ac);
        }
        fprintf(m_pLog, "%-30s\t%-14s\t%-14s\t%s\n","...","...","...","...");
        // last 5
        for(int i = (n-5); i < n; i++)
        {
            fprintf(m_pLog, "%-30s\t%f\t%f\t%f\n", data[i].name, data[i].x, data[i].y, data[i].area_ac);
        }
    }
    fprintf(m_pLog, "\n[%d rows x 3 columns]\n", n);
    fprintf(m_pLog, "\n");
}/* end PrintSubcatchmentData() */

/* ----------------------------------------------------------------------------
   PrintClusteredData()

    A function that displays subcatchment data along with clustering results.
---------------------------------------------------------------------------- */
void KMeansParamInitializer::PrintClusteredData(SubcatchmentCoordStruct * data, int n)
{
    fprintf(m_pLog, "%-5s\t%-30s\t%-14s\t%-14s\t%-9s\t%-7s\n",
           "index", "NAME", "CX", "CY", "subc_Area", "cluster");
    if(n < 10)
    {
        for(int i = 0; i < n; i++)
        {
            fprintf(m_pLog, "%-5d\t%-30s\t%f\t%f\t%f\t%7d\n", 
            (i+1), data[i].name, data[i].x, data[i].y, data[i].area_ac, data[i].cluster_id);
        }
    }
    else
    {
        // first 5
        for(int i = 0; i < 5; i++)
        {
            fprintf(m_pLog, "%-5d\t%-30s\t%f\t%f\t%f\t%7d\n", 
            (i+1), data[i].name, data[i].x, data[i].y, data[i].area_ac, data[i].cluster_id);
        }
        fprintf(m_pLog, "%-5s\t%-30s\t%-14s\t%-14s\t%-9s\t%7s\n",
               "...","...","...","...","...","...");
        // last 5
        for(int i = (n-5); i < n; i++)
        {
            fprintf(m_pLog, "%-5d\t%-30s\t%f\t%f\t%f\t%7d\n", 
            (i+1), data[i].name, data[i].x, data[i].y, data[i].area_ac, data[i].cluster_id);
        }
    }
    fprintf(m_pLog, "\n[%d rows x 5 columns]\n", n);
    fprintf(m_pLog, "\n");
}/* end PrintClusteredData() */

/* ----------------------------------------------------------------------------
   PrintWithAreas())

    A function that displays subcatchment data along with cluster areas.
---------------------------------------------------------------------------- */
void KMeansParamInitializer::PrintWithAreas(SubcatchmentCoordStruct * data, int n)
{
    fprintf(m_pLog, "%-5s\t%-30s\t%-14s\t%-14s\t%-9s\t%-7s\t%-12s\n",
           "index", "NAME", "CX", "CY", "subc_Area", "cluster", "Cluster_area");
    if(n < 10)
    {
        for(int i = 0; i < n; i++)
        {
            fprintf(m_pLog, "%-5d\t%-30s\t%f\t%f\t%f\t%7d\t%12.6f\n", 
            (i+1), data[i].name, data[i].x, data[i].y, data[i].area_ac, data[i].cluster_id, data[i].cluster_area_ac);
        }
    }
    else
    {
        // first 5
        for(int i = 0; i < 5; i++)
        {
            fprintf(m_pLog, "%-5d\t%-30s\t%f\t%f\t%f\t%7d\t%12.6f\n", 
            (i+1), data[i].name, data[i].x, data[i].y, data[i].area_ac, data[i].cluster_id, data[i].cluster_area_ac);
        }
        fprintf(m_pLog, "%-5s\t%-30s\t%-14s\t%-14s\t%-9s\t%7s\t%12s\n",
               "...","...","...","...","...","...","...");
        // last 5
        for(int i = (n-5); i < n; i++)
        {
            fprintf(m_pLog, "%-5d\t%-30s\t%f\t%f\t%f\t%7d\t%12.6f\n", 
            (i+1), data[i].name, data[i].x, data[i].y, data[i].area_ac, data[i].cluster_id, data[i].cluster_area_ac);
        }
    }
    fprintf(m_pLog, "\n[%d rows x 6 columns]\n", n);
    fprintf(m_pLog, "\n");
} /* end PrintWithAreas() */

/* ----------------------------------------------------------------------------
   ReadCoordFile()

    A function that reads in a subcatfchment file where columns are separated 
    by 'sep'
---------------------------------------------------------------------------- */
SubcatchmentCoordStruct * KMeansParamInitializer::ReadCoordsFile(FILE * pIn, char sep, int * n)
{
    char line[1000];
    SubcatchmentCoordStruct * pSubInfo;
    int i, j;
    char * token;
    char s[2] = {sep, 0x0};
    char * header = "NAME	X	Y	Area";
    *n = 0;

    rewind(pIn);
    fgets(line, 1000, pIn);
    rtrim(line);
    // validate header
    if(strcmp(line, header) != 0)
    {
        printf("Error - coords file lacks header!\n");
        printf("Expected: |%s|\n", header);
        printf("Actual  : |%s|\n", line);
        return NULL;
    }
    // count number of entries
    while(1)
    {
        fgets(line, 1000, pIn);
        if(! feof(pIn))
        {
            (*n)++;
        }
        else
        {
            break;
        }
    }/* end while() */
    if((*n) > 0)
    {
        pSubInfo = new SubcatchmentCoordStruct[(*n)];
        rewind(pIn);
        fgets(line, 1000, pIn);
        for(i = 0; i < (*n); i++)
        {
            fgets(line, 1000, pIn);
            token = strtok(line, s);
            j = 0;
            while(token != NULL)
            {
                if(j == 0)
                {
                    strncpy(pSubInfo[i].name, token, 1000);
                }
                else if(j == 1)
                {
                    pSubInfo[i].x = atof(token);
                }
                else if(j == 2)
                {
                    pSubInfo[i].y = atof(token);
                }
                else if(j == 3)
                {
                    pSubInfo[i].area_ac = atof(token);
                }
                j++;
                token = strtok(NULL, s);
            }/* end while() */
            if(j != 4)
            {
                printf("Parsing error in row %d\n", i+1);
                printf("Expected 4 fields and got %d fields instead.\n", j);
                (*n) = i + 1;
                break;
            }
        }/* end for() */
    }
    else
    {
        pSubInfo = NULL;
       *n = 0;
    }
    
    return pSubInfo;
} /* end ReadCoordsFile() */

/******************************************************************************
GetParameterSets()

Apply the algorithm and populate the "pVals" matrix with initial parameter sets.
******************************************************************************/
void KMeansParamInitializer::GetParameterSets(double ** pVals, int start)
{
   int i, j, end;

   end = start + m_NumSets;
   for(i = start; i < end; i++)
   {
      for(j = 0; j < m_NumParams; j++)
      {
         pVals[i][j] = m_pParams->GetParamPtr(j)->ConvertInVal(atof("0.0"));
      }
   }
	return; 
} /* end GetParameterSets() */

/******************************************************************************
CTOR
  
Assign member variables.
******************************************************************************/
KMeansParamInitializer::KMeansParamInitializer(ParameterGroup * pParamGroup, FILE * pInFile)
{
   int i, j, k, num;
   char * pTok;
   char * line;
   char tmp[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];
   char pFileName[] = "ostIn.txt";

   m_pScenarios = NULL;
   m_pScenarioClusters = NULL;    
   m_pDataOriginal = NULL;    
   m_KMEANS.labels = NULL;
   m_KMEANS.centers = NULL;

   m_NumSets = 10;
   m_NumClusters = 0;
   m_pParams = pParamGroup;
   m_NumParams = pParamGroup->GetNumParams();
   strcpy(m_LogFileName, "OstKMeansParamInit.log");
   strcpy(m_CoordsFile, "CSO014_ParkingSAimp_coord_andAREA.txt");
   m_AreaPerScenario = 50.0;
   m_AreaPerLID = 440.0;
   m_AreaConversionFactor = 43560.0;

   m_pLog = fopen(m_LogFileName, "w");

   //read in KMeans configuration
   if(pInFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open KMeans config. file. Using Defaults");
      fclose(m_pLog);      
      return;
   }/* end if() */

   rewind(pInFile);

   //make sure correct tokens are present
   if(CheckToken(pInFile, "BeginKMeansInitializer", pFileName) == true)
   {
      FindToken(pInFile, "EndKMeansInitializer", pFileName);
      rewind(pInFile);

      FindToken(pInFile, "BeginKMeansInitializer", pFileName);
      line = GetNxtDataLine(pInFile, pFileName);
      while(strstr(line, "EndKMeansInitializer") == NULL)
      {         
         if(strstr(line, "AreaPerScenario") != NULL)
         {
            sscanf(line, "%s %f", tmp, &m_AreaPerScenario); 
         }/*end else if() */
         else if(strstr(line, "AreaPerLID") != NULL)
         {
            sscanf(line, "%s %f", tmp, &m_AreaPerLID); 
         }/*end else if() */         
         else if(strstr(line, "AreaConversionFactor") != NULL)
         {
            sscanf(line, "%s %f", tmp, &m_AreaConversionFactor); 
         }/*end else if() */         
         else if(strstr(line, "CoordsFile") != NULL)
         {
            sscanf(line, "%s %s", tmp, m_CoordsFile); 
         }/*end else if() */
         else if(strstr(line, "NumSets") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumSets); 
         }/*end else if() */
         else if(strstr(line, "NumScenarios") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumSets); 
         }/*end else if() */
         else if(strstr(line, "NumClusters") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumClusters); 
         }/*end else if() */         
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pInFile, pFileName);
      } /* end while() */
   }/* end if() */   

   rewind(pInFile);

   Configure();

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()

Free up memory of the initializer.
******************************************************************************/
void KMeansParamInitializer::Destroy(void)
{
    // free up dynamic memory
    if(m_pScenarios != NULL)
    { 
      for(int i = 0; i < m_NumSets; i++)
      {
         delete [] m_pScenarios[i];
      }
      delete [] m_pScenarios;
    }
    if(m_pScenarioClusters != NULL)
    { 
      for(int i = 0; i < m_NumSets; i++)
      {
         delete [] m_pScenarioClusters[i];
      }
      delete [] m_pScenarioClusters;
    }
    if(m_pDataOriginal != NULL)
    {
      delete [] m_pDataOriginal;
    }
    if(m_KMEANS.labels != NULL)
    {
      delete [] m_KMEANS.labels;
    }
    if(m_KMEANS.centers != NULL)
    {
      delete [] m_KMEANS.centers;
    }
   IncDtorCount();
   fclose(m_pLog);
} /* end Destroy() */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void KMeansParamInitializer::Write(FILE * pFile, int type)
{
   fprintf(pFile, "********** Parameter Initialization **********\n");
   fprintf(pFile, "Name                   : KMeans Clustering\n");
   fprintf(pFile, "Log File               : %s\n", m_LogFileName);
   fprintf(pFile, "Coordinates File       : %s\n", m_CoordsFile);
   fprintf(pFile, "Num Scenarios          : %d\n", m_NumSets);
   fprintf(pFile, "Area per Scenario      : %f\n", m_AreaPerScenario);
   fprintf(pFile, "Area per LID           : %f\n", m_AreaPerLID);
   fprintf(pFile, "Area Conversion Factor : %f\n", m_AreaConversionFactor);
   fprintf(pFile, "Num Clusters           : %d\n", m_NumClusters);
} /* end Write() */

/******************************************************************************
Configure()

Read coordinate data from input file, setup KMeans clusters, and map to 
parameter sets.
******************************************************************************/
int KMeansParamInitializer::Configure(void)
{
    SubcatchmentCoordStruct * data_w_clusters; // data with cluster id added
    SubcatchmentCoordStruct * data_merged; // data with cluster id and cluster areas added

    fprintf(m_pLog, "# Using k-means clustering for LIDs allocation\n\n");
    fprintf(m_pLog,"## Introduction\n");
    fprintf(m_pLog,"Previous research states that aggregating LIDs is a good strategy \n");
    fprintf(m_pLog,"to reduce runoff and Combined Sewer Overflows.\n\n");
    fprintf(m_pLog,"In this document we will exploite the k-means clustering algorithm \n");
    fprintf(m_pLog,"to test wether aggregation increase the performance of Pervious \n");
    fprintf(m_pLog,"Pavements to reduce CSO.\n\n");

    fprintf(m_pLog,"### Step 1: Cluster the subcatchments by location\n");
    fprintf(m_pLog,"We need the set of coordinates for each subcatchment\n\n");

    // Create clusters from LID coordinates
    int clusters_set = m_NumClusters; // number of desired clusters

    // Load the text file contains the xy coordinates and areas
    FILE * pCoordFile = fopen(m_CoordsFile, "r");
    if(pCoordFile == NULL)
    {
        fprintf(m_pLog, "Error - can't open coordinate file (CSO014_ParkingSAimp_coord_andAREA.txt)\n");
        return -1;
    }
    // read the coordinate file
    int n_subcatchments = 0;
    m_pDataOriginal = ReadCoordsFile(pCoordFile,'\t', &n_subcatchments);
    fclose(pCoordFile);

    // assign the index to be the subcatchment name
    char data_original_index[] = "NAME";

    // subcatchm ent area is given in ac
    PrintSubcatchmentData(m_pDataOriginal, n_subcatchments);

    // Create clusters and add a column to the original dataset
    m_KMEANS.num_clusters = clusters_set;
    m_KMEANS.num_data_points = n_subcatchments;
    m_KMEANS.labels = new int[m_KMEANS.num_data_points];
    m_KMEANS.centers = new KMeansCoordStruct[m_KMEANS.num_clusters];
    CreateClusters(m_pDataOriginal, &m_KMEANS); 
    KMeansCoordStruct * centers = m_KMEANS.centers;
    int * labels = m_KMEANS.labels; // which cluster each coord belongs to
    data_w_clusters = AppendLabels(m_pDataOriginal , labels, n_subcatchments);
    PrintClusteredData(data_w_clusters, n_subcatchments);

    // Calculate the area assigned to each cluster
    data_merged = AppendClusterAreas(data_w_clusters, n_subcatchments, m_KMEANS.num_clusters);
    PrintWithAreas(data_merged, n_subcatchments);

    fprintf(m_pLog, "### Step 2: Create scenarios for LID location based on clusters\n\n");
    int i;
    int j;
    int r;
    int s;
    int r_count;
    int s_count;
    int num_scenarios = m_NumSets;
    double area_per_scenario = m_AreaPerScenario;
    double area_per_lid = m_AreaPerLID;
    double area_conversion_factor = m_AreaConversionFactor;
    double assigned_area;
    double cluster_area;
    m_pScenarios = new LidScenarioStruct *[num_scenarios];
    m_pScenarioClusters = new int *[num_scenarios];
    for(i = 0; i < num_scenarios; i++)
    {
        m_pScenarios[i] = new LidScenarioStruct[n_subcatchments];
        m_pScenarioClusters[i] = new int[m_KMEANS.num_clusters];
        for(j = 0; j < n_subcatchments; j++)
        {
            m_pScenarios[i][j].num_lids = 0;
            m_pScenarios[i][j].name = data_merged[j].name;
        }
        for(j = 0; j < m_KMEANS.num_clusters; j++)
        {
            m_pScenarioClusters[i][j] = -1;
        }
    }
    r_count = 0;
    for(i = 0; i < num_scenarios; i++)
    {
        assigned_area = 0;
        s_count = 0;
        while(assigned_area < area_per_scenario)
        {
            // pick a random cluster and determine how much area can be assigned to it
            r_count++;
            r = PickRandomCluster(r_count);
            m_pScenarioClusters[i][s_count] = r;
            s_count++;
            cluster_area = GetClusterArea(data_merged, n_subcatchments, r);
            if((assigned_area + cluster_area) > area_per_scenario)
            {
               cluster_area = (area_per_scenario - assigned_area);
            }
            // update amount of assigned area
            assigned_area += cluster_area;

            // assign LIDs to the r-th cluster until cluster_area is used up ....
            for(j = 0; j < n_subcatchments; j++)
            {
                if(data_merged[j].cluster_id == r)
                {
                    if(data_merged[j].area_ac > cluster_area)
                    {
                        m_pScenarios[i][j].num_lids = ((cluster_area * area_conversion_factor) / area_per_lid);
                        cluster_area = 0.00;
                    }
                    else
                    {
                        m_pScenarios[i][j].num_lids = (((data_merged[j].area_ac) * area_conversion_factor) / area_per_lid);
                        cluster_area = cluster_area - data_merged[j].area_ac;
                    }
                }/* end if(member of r-th cluster) */
                if(cluster_area <= 0.01)
                {
                    break;
                }
            }/* end for(each sub_catchment) */
        }/* end while(assigned_area) */
    }/* end for(each scenario) */

    // display results of lid assignment
    fprintf(m_pLog, "%-4s\t%-30s\t%-7s\t%-8s\n", "iter", "Subcatchment_Name", "Cluster", "num_LIDs");
    for(i = 0; i < num_scenarios; i++)
    {
        for(s = 0; s < m_KMEANS.num_clusters; s++)
        {
            r = m_pScenarioClusters[i][s];
            if(r != -1)
            {
                for(j = 0; j < n_subcatchments; j++)
                {
                    if((m_pScenarios[i][j].num_lids > 0) && (data_merged[j].cluster_id == r))
                    {
                        fprintf(m_pLog, "%-4d\t%-30s\t%-7d\t%8.2f\n", i, m_pScenarios[i][j].name,  
                                           data_merged[j].cluster_id, m_pScenarios[i][j].num_lids);
                    }
                }
                fprintf(m_pLog, "\n");
            }
        }
    }
    
    return 0;
}/* end Configure()) */
