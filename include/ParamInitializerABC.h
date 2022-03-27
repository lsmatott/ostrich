/******************************************************************************
File      : ParamInitializerABC.h
Author    : L. Shawn Matott 
Copyright : 2022, L. Shawn Matott

This is an abstract base class that supports methods for assigning sets of
initial parameter values.

Version History
02-27-2022   lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef PARAM_INITIALIZER_ABC_H
#define PARAM_INITIALIZER_ABC_H

#include "MyHeaderInc.h"

//forward decs
class ParameterGroup;
class ParameterABC;

#define PMAP_EXPLICIT (0)
#define PMAP_IMPLICIT (1)

/******************************************************************************
class ParamInitializerABC (parameter initializer) base class
******************************************************************************/
class ParamInitializerABC
{
    public :
    	virtual ~ParamInitializerABC(void){ DBG_PRINT("ParamInitializerABC::DTOR"); }
		virtual void Destroy(void)=0;
		virtual void GetParameterSets(double ** pVals, int start)=0;
		virtual void Write(FILE * pFile, int type)=0; 
		virtual int GetNumParameterSets(void)=0;
}; /* end class ParamInitializerABC */

/******************************************************************************
class HamedParamInitializer

Initialize sets of parameter values using Hamed's method. This is thought to be 
useful for capturing larger sections of the pareto front in multi-objective problems.
******************************************************************************/
class HamedParamInitializer : public ParamInitializerABC
{
    public:
	~HamedParamInitializer(void){ DBG_PRINT("HamedParamInitializer::DTOR"); Destroy(); }
	void Destroy(void);
	HamedParamInitializer(ParameterGroup * pParamGroup, FILE * pInFile);
	void GetParameterSets(double ** pVals, int start);
	int GetNumParameterSets(void){ return m_NumSets; }
	void Write(FILE * pFile, int type);

private:
    void GetSample(double * stest,  int sindex);
	ParameterGroup * m_pParams;
	int m_HamedOffset;
	int m_NumSets;
	int m_NumParams;
}; /* end class HamedParamInitializer */

/******************************************************************************
Supporting data structures for KMeansParamInitializer
******************************************************************************/
typedef struct SUBCATCHMENT_COORD_STRUCT {
   char name[1000];
   char pname[1000];
   double x;
   double y;
   double area_ac;
   int cluster_id;
   double cluster_area_ac;
} SubcatchmentCoordStruct;

typedef struct KMEANS_COORD_STRUCT {
    double x;
    double y;
}KMeansCoordStruct;

typedef struct KMEANS_CLUSTER_STRUCT {
   int num_clusters;
   int num_data_points; 
   int * labels; // which cluster each data point belongs to
   KMeansCoordStruct * centers; // the center point of each cluster
}KMeansClusterStruct;

typedef struct LID_SCENARIO {
   char * name;
   char * pname;
   double num_lids;
} LidScenarioStruct;

/******************************************************************************
class KMeansParamInitializer

Initialize parameter values using K-Means clustering. This is thought to be useful
for grouping spatially related parameters.
******************************************************************************/
class KMeansParamInitializer : public ParamInitializerABC
{
    public:
	~KMeansParamInitializer(void){ DBG_PRINT("KMeansParamInitializer::DTOR"); Destroy(); }
	void Destroy(void);
	KMeansParamInitializer(ParameterGroup * pParamGroup, FILE * pInFIle);
	void GetParameterSets(double ** pVals, int start);
	void Write(FILE * pFile, int type);
	int GetNumParameterSets(void){ return m_NumSets; }

private:
    int Configure(void);
    double GetClusterArea(SubcatchmentCoordStruct * pData, int n, int cluster_id);
    int PickRandomCluster(int n);
    SubcatchmentCoordStruct * AppendClusterAreas(SubcatchmentCoordStruct * data, int n_subcatchments, int n_clusters);
    SubcatchmentCoordStruct * AppendLabels(SubcatchmentCoordStruct * data, int * labels, int n);
    void CreateClusters(SubcatchmentCoordStruct * data, KMeansClusterStruct * pKMeans);
	SubcatchmentCoordStruct * ReadCoordsFile(FILE * pIn, char sep, int * n);
	void PrintWithAreas(SubcatchmentCoordStruct * data, int n);
    void PrintClusteredData(SubcatchmentCoordStruct * data, int n);
    void PrintSubcatchmentData(SubcatchmentCoordStruct * data, int n);
    void rtrim(char * pStr);
    void AssignImplicitParameterNames(void);

    KMeansClusterStruct m_KMEANS;
    SubcatchmentCoordStruct * m_pDataOriginal;
    LidScenarioStruct ** m_pScenarios;
    int ** m_pScenarioClusters;

	ParameterGroup * m_pParams;
    int m_ParameterMapping; /* how to map subcatchment names to parameter names */
    int m_NumSets; /* num scenarios */
	int m_NumParams;
	int m_NumClusters;
    double m_AreaPerScenario;
    double m_AreaPerLID;
    double m_AreaConversionFactor;
	char m_LogFileName[DEF_STR_SZ];
	char m_CoordsFile[DEF_STR_SZ];
	FILE * m_pLog;
}; /* end class KMeansParamInitializer */

#endif /* PARAM_INITIALIZER_ABC_H */






