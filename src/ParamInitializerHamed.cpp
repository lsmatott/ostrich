/******************************************************************************
File      : HamedParamInitializer.cpp
Author    : L. Shawn Matott
Copyright : 2022, L. Shawn Matott

Encapsulates the Hamed parameter initialization strategy.

Version History
02-27-2022    lsm   created
******************************************************************************/
#include <string.h>
#include "ParamInitializerABC.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "Exception.h"
#include "Utility.h"
#include "mpi.h"

/******************************************************************************
GetParameterSets()

Apply the algorithm and populate the "pVals" matrix with initial parameter values.
******************************************************************************/
void HamedParamInitializer::GetParameterSets(double ** pVals, int start)
{
   int i, j, end;

   end = start + m_NumSets;
   for(i = start; i < end; i++)
   {
      GetSample(pVals[i], (i - start));
   }
	return; 
} /* end GetParameterSets() */

/******************************************************************************
CTOR
  
Assign member variables by parsing input file.
******************************************************************************/
HamedParamInitializer::HamedParamInitializer(ParameterGroup * pParamGroup, FILE * pInFile)
{
   int i, j, k, num;
   char * pTok;
   char * line;
   char tmp[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];
   char pFileName[] = "ostIn.txt";

   m_HamedOffset = 397;
   MPI_Comm_size(MPI_COMM_WORLD, &m_NumSets);
   m_pParams = pParamGroup;
   m_NumParams = pParamGroup->GetNumParams();

   //read in Hamed configuration
   if(pInFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open Hamed config. file. Using Defaults");      
      return;
   }/* end if() */

   rewind(pInFile);

   //make sure correct tokens are present
   if(CheckToken(pInFile, "BeginHamedInitializer", pFileName) == true)
   {
      FindToken(pInFile, "EndHamedInitializer", pFileName);
      rewind(pInFile);

      FindToken(pInFile, "BeginHamedInitializer", pFileName);
      line = GetNxtDataLine(pInFile, pFileName);
      while(strstr(line, "EndHamedInitializer") == NULL)
      {         
         if(strstr(line, "HamedOffset") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_HamedOffset); 
         }/*end else if() */         
         else if(strstr(line, "NumSets") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumSets); 
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

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
InitUsingHamedMethod()

Initialize the solutions using Hamed's carefully constructed empirical method.
******************************************************************************/
void HamedParamInitializer::GetSample(double * stest,  int sindex)
{
   int n, offset, nsubs;
   double pminj, pmaxj;
   double pminjo, pmaxjo;

   nsubs = m_NumSets;
   n = m_NumParams;
   offset = m_HamedOffset;

   // first pass - use urandom
   for (int j = 0; j < n; j++)
   {
      pminj = m_pParams->GetParamPtr(j)->GetLwrBnd();
      pmaxj = m_pParams->GetParamPtr(j)->GetUprBnd();
      stest[j] = pminj + (pmaxj - pminj)*UniformRandom();
   }/* end for() */

   // second pass - update selected entries using Hamed's empirical method
   for (int j = 0; j < n/3; j++)
   {
      pminj = m_pParams->GetParamPtr(j)->GetLwrBnd();
      pmaxj = m_pParams->GetParamPtr(j)->GetUprBnd();
      pminjo = m_pParams->GetParamPtr(j+offset)->GetLwrBnd();
      pmaxjo = m_pParams->GetParamPtr(j+offset)->GetUprBnd();

      // The if structure distributes the initial conditions across the range of the "offset" term
      if (((double)sindex/(double)nsubs) <= 0.1) 
      {                           
         stest[j] = pminj;
         if ( (((double)sindex/(double)nsubs) >= 0.0) && (((double)sindex/(double)nsubs) <= 0.01) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.01) && (((double)sindex/(double)nsubs) <= 0.02)) {
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.02) && (((double)sindex/(double)nsubs) <= 0.03)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.03) && (((double)sindex/(double)nsubs) <= 0.04)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.04) && (((double)sindex/(double)nsubs) <= 0.05)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.05) && (((double)sindex/(double)nsubs) <= 0.06)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.06) && (((double)sindex/(double)nsubs) <= 0.07)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.07) && (((double)sindex/(double)nsubs) <= 0.08)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.08) && (((double)sindex/(double)nsubs) <= 0.09)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.09) && (((double)sindex/(double)nsubs) <= 0.098)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.098) && (((double)sindex/(double)nsubs) <= 0.10)){
            stest[j+offset] = pmaxjo; }
     } 
     else if ( (((double)sindex/(double)nsubs) > 0.1) && (((double)sindex/(double)nsubs) <= 0.2) ) 
     {
         stest[j] = pminj + (pmaxj - pminj) * 0.1;
         if ( (((double)sindex/(double)nsubs) > 0.10) && (((double)sindex/(double)nsubs) <= 0.11) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.11) && (((double)sindex/(double)nsubs) <= 0.12)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.12) && (((double)sindex/(double)nsubs) <= 0.13)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.13) && (((double)sindex/(double)nsubs) <= 0.14)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.14) && (((double)sindex/(double)nsubs) <= 0.15)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.15) && (((double)sindex/(double)nsubs) <= 0.16)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.16) && (((double)sindex/(double)nsubs) <= 0.17)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.17) && (((double)sindex/(double)nsubs) <= 0.18)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.18) && (((double)sindex/(double)nsubs) <= 0.19)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.19) && (((double)sindex/(double)nsubs) <= 0.198)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.198) && (((double)sindex/(double)nsubs) <= 0.20)){
            stest[j+offset] = pmaxjo; }       
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.2) && (((double)sindex/(double)nsubs) <= 0.3) ) 
      {
         stest[j] = pminj + (pmaxj - pminj) * 0.2;
         if ( (((double)sindex/(double)nsubs) > 0.20) && (((double)sindex/(double)nsubs) <= 0.21) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.21) && (((double)sindex/(double)nsubs) <= 0.22)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.22) && (((double)sindex/(double)nsubs) <= 0.23)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.23) && (((double)sindex/(double)nsubs) <= 0.24)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.24) && (((double)sindex/(double)nsubs) <= 0.25)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.25) && (((double)sindex/(double)nsubs) <= 0.26)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.26) && (((double)sindex/(double)nsubs) <= 0.27)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.27) && (((double)sindex/(double)nsubs) <= 0.28)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.28) && (((double)sindex/(double)nsubs) <= 0.29)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.29) && (((double)sindex/(double)nsubs) <= 0.298)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.298) && (((double)sindex/(double)nsubs) <= 0.30)){
            stest[j+offset] = pmaxjo; }
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.3) && (((double)sindex/(double)nsubs) <= 0.4) ) 
      {
         stest[j] = pminj + (pmaxj - pminj) * 0.3; 
         if ( (((double)sindex/(double)nsubs) > 0.30) && (((double)sindex/(double)nsubs) <= 0.31) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.31) && (((double)sindex/(double)nsubs) <= 0.32)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.32) && (((double)sindex/(double)nsubs) <= 0.33)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.33) && (((double)sindex/(double)nsubs) <= 0.34)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.34) && (((double)sindex/(double)nsubs) <= 0.35)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.35) && (((double)sindex/(double)nsubs) <= 0.36)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.36) && (((double)sindex/(double)nsubs) <= 0.37)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.37) && (((double)sindex/(double)nsubs) <= 0.38)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.38) && (((double)sindex/(double)nsubs) <= 0.39)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.39) && (((double)sindex/(double)nsubs) <= 0.398)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.398) && (((double)sindex/(double)nsubs) <= 0.40)){
            stest[j+offset] = pmaxjo; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.4) && (((double)sindex/(double)nsubs) <= 0.5) ) 
      {
         stest[j] = pminj + (pmaxj - pminj) * 0.4;
         if ( (((double)sindex/(double)nsubs) > 0.40) && (((double)sindex/(double)nsubs) <= 0.41) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.41) && (((double)sindex/(double)nsubs) <= 0.42)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.42) && (((double)sindex/(double)nsubs) <= 0.43)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.43) && (((double)sindex/(double)nsubs) <= 0.44)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.44) && (((double)sindex/(double)nsubs) <= 0.45)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.45) && (((double)sindex/(double)nsubs) <= 0.46)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.46) && (((double)sindex/(double)nsubs) <= 0.47)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.47) && (((double)sindex/(double)nsubs) <= 0.48)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.48) && (((double)sindex/(double)nsubs) <= 0.49)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.49) && (((double)sindex/(double)nsubs) <= 0.498)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.498) && (((double)sindex/(double)nsubs) <= 0.50)){
            stest[j+offset] = pmaxjo; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.5) && (((double)sindex/(double)nsubs) <= 0.6) ) 
      {
         stest[j] = pminj + (pmaxj - pminj) * 0.5;
         if ( (((double)sindex/(double)nsubs) > 0.50) && (((double)sindex/(double)nsubs) <= 0.51) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.51) && (((double)sindex/(double)nsubs) <= 0.52)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.52) && (((double)sindex/(double)nsubs) <= 0.53)){
             stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.53) && (((double)sindex/(double)nsubs) <= 0.54)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.54) && (((double)sindex/(double)nsubs) <= 0.55)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.55) && (((double)sindex/(double)nsubs) <= 0.56)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.56) && (((double)sindex/(double)nsubs) <= 0.57)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.57) && (((double)sindex/(double)nsubs) <= 0.58)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.58) && (((double)sindex/(double)nsubs) <= 0.59)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.59) && (((double)sindex/(double)nsubs) <= 0.598)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.598) && (((double)sindex/(double)nsubs) <= 0.60)){
            stest[j+offset] = pmaxjo; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.6) && (((double)sindex/(double)nsubs) <= 0.7) ) 
      {
         stest[j] = pminj + (pmaxj - pminj) * 0.6;
         if ( (((double)sindex/(double)nsubs) > 0.60) && (((double)sindex/(double)nsubs) <= 0.61) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.61) && (((double)sindex/(double)nsubs) <= 0.62)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.62) && (((double)sindex/(double)nsubs) <= 0.63)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.63) && (((double)sindex/(double)nsubs) <= 0.64)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.64) && (((double)sindex/(double)nsubs) <= 0.65)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.65) && (((double)sindex/(double)nsubs) <= 0.66)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.66) && (((double)sindex/(double)nsubs) <= 0.67)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.67) && (((double)sindex/(double)nsubs) <= 0.68)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.68) && (((double)sindex/(double)nsubs) <= 0.69)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.69) && (((double)sindex/(double)nsubs) <= 0.698)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.698) && (((double)sindex/(double)nsubs) <= 0.70)){
            stest[j+offset] = pmaxjo; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.7) && (((double)sindex/(double)nsubs) <= 0.8) ) 
      {
         stest[j] = pminj + (pmaxj - pminj) * 0.7;
         if ( (((double)sindex/(double)nsubs) > 0.70) && (((double)sindex/(double)nsubs) <= 0.71) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.71) && (((double)sindex/(double)nsubs) <= 0.72)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.72) && (((double)sindex/(double)nsubs) <= 0.73)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.73) && (((double)sindex/(double)nsubs) <= 0.74)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.74) && (((double)sindex/(double)nsubs) <= 0.75)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.75) && (((double)sindex/(double)nsubs) <= 0.76)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.76) && (((double)sindex/(double)nsubs) <= 0.77)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.77) && (((double)sindex/(double)nsubs) <= 0.78)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.78) && (((double)sindex/(double)nsubs) <= 0.79)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.79) && (((double)sindex/(double)nsubs) <= 0.798)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.798) && (((double)sindex/(double)nsubs) <= 0.80)){
            stest[j+offset] = pmaxjo; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.8) && (((double)sindex/(double)nsubs) <= 0.9) ) 
      {
         stest[j] = pminj + (pmaxj - pminj) * 0.8;
         if ( (((double)sindex/(double)nsubs) > 0.80) && (((double)sindex/(double)nsubs) <= 0.81) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.81) && (((double)sindex/(double)nsubs) <= 0.82)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.1; 
         }else if((((double)sindex/(double)nsubs) > 0.82) && (((double)sindex/(double)nsubs) <= 0.83)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.2; 
         }else if((((double)sindex/(double)nsubs) > 0.83) && (((double)sindex/(double)nsubs) <= 0.84)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.3; 
         }else if((((double)sindex/(double)nsubs) > 0.84) && (((double)sindex/(double)nsubs) <= 0.85)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.4; 
         }else if((((double)sindex/(double)nsubs) > 0.85) && (((double)sindex/(double)nsubs) <= 0.86)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.5; 
         }else if((((double)sindex/(double)nsubs) > 0.86) && (((double)sindex/(double)nsubs) <= 0.87)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.6; 
         }else if((((double)sindex/(double)nsubs) > 0.87) && (((double)sindex/(double)nsubs) <= 0.88)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.7; 
         }else if((((double)sindex/(double)nsubs) > 0.88) && (((double)sindex/(double)nsubs) <= 0.89)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.8; 
         }else if((((double)sindex/(double)nsubs) > 0.89) && (((double)sindex/(double)nsubs) <= 0.898)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.9; 
         }else if((((double)sindex/(double)nsubs) > 0.898) && (((double)sindex/(double)nsubs) <= 0.90)){
            stest[j+offset] = pmaxjo; }   
      } 
      else if ( (((double)sindex/(double)nsubs) > 0.9) && (((double)sindex/(double)nsubs) <= 0.95) ) 
      {
         stest[j] = pminj + (pmaxj - pminj) * 0.9;
         if ( (((double)sindex/(double)nsubs) > 0.90) && (((double)sindex/(double)nsubs) <= 0.91) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.91) && (((double)sindex/(double)nsubs) <= 0.92)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.25; 
         }else if((((double)sindex/(double)nsubs) > 0.92) && (((double)sindex/(double)nsubs) <= 0.93)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.50; 
         }else if((((double)sindex/(double)nsubs) > 0.93) && (((double)sindex/(double)nsubs) <= 0.94)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.75; 
         }else if((((double)sindex/(double)nsubs) > 0.94) && (((double)sindex/(double)nsubs) <= 0.95)){
            stest[j+offset] = pmaxjo; }   
      } 
      else 
      {
         stest[j] = pmaxj;
         if ( (((double)sindex/(double)nsubs) > 0.95) && (((double)sindex/(double)nsubs) <= 0.96) ) {
            stest[j+offset] = pminjo; 
         }else if((((double)sindex/(double)nsubs) > 0.96) && (((double)sindex/(double)nsubs) <= 0.97)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.25; 
         }else if((((double)sindex/(double)nsubs) > 0.97) && (((double)sindex/(double)nsubs) <= 0.98)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.50; 
         }else if((((double)sindex/(double)nsubs) > 0.98) && (((double)sindex/(double)nsubs) <= 0.99)){
            stest[j+offset] = pminjo + (pmaxjo - pminjo) * 0.75; 
         }else if((((double)sindex/(double)nsubs) > 0.99) && (((double)sindex/(double)nsubs) <= 1)){
            stest[j+offset] = pmaxjo; }
      }       
   } /* end for() */
} /* end InitUsingHamedMethod() */

/******************************************************************************
Destroy()

Free up memory of the initializer.
******************************************************************************/
void HamedParamInitializer::Destroy(void)
{
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
Write()

Writes formatted output to pFile.
******************************************************************************/
void HamedParamInitializer::Write(FILE * pFile, int type)
{
   fprintf(pFile, "****** Parameter Initialization ******\n");
   fprintf(pFile, "Name       : Hamed's Method\n");
   fprintf(pFile, "Offset     : %d\n", m_HamedOffset);
   fprintf(pFile, "Num Sets   : %d\n", m_NumSets);
} /* end Write() */

