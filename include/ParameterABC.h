/******************************************************************************
File      : ParameterABC.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a parameter. Parameters are variables in the model which are to 
be calibrated.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
11-25-03    lsm   Modified to support multiple stages of unit conversion.
07-05-04    lsm   added integer and combinatorial parameter support
12-02-04    lsm   made ConvertInVal() a required member function
03-03-05    lsm   Added support for ON/OFF parameter threshold
******************************************************************************/
#ifndef PARAMETER_ABC_H
#define PARAMETER_ABC_H

#include "MyHeaderInc.h"
#include "string"

// forward decs
class ConstraintABC;

/* 
An enum is defined for each type of transformation.
*/
typedef enum TRANSFORM_TYPE
{
   TX_NONE   = 0,
   TX_LOG10  = 1,
   TX_LN     = 2
}TransformTypeEnum;

/* 
An enum is defined for each stage of transformation.
*/
typedef enum TRANSFORM_STAGE
{
   TX_IN   = 0,
   TX_OST  = 1,
   TX_OUT  = 2
}TransformStageEnum;
#define NUM_STAGES (3)

/******************************************************************************
class ParameterABC

Abstract base class of a parameter.
******************************************************************************/
class ParameterABC
{
   public:       
      virtual ~ParameterABC(void){ DBG_PRINT("ParameterABC::DTOR"); }     
      virtual void Destroy(void) = 0;
      virtual void   GetValAsStr(UnmoveableString valStr) = 0;
      virtual void   Write(FILE * pFile, int type) = 0;
      virtual double GetInitialValueTransformed(void) = 0;
      virtual double GetLowerBoundTransformed(void) = 0;
      virtual double GetUpperBoundTransformed(void) = 0;
      virtual void SetLowerBoundTransformed(double val) = 0;
      virtual void SetUpperBoundTransformed(double val) = 0;
      virtual double GetEstimatedValueTransformed(void) = 0;
      virtual double SetEstimatedValueTransformed(double estVal) = 0;

      //threshold values (allow for implicit on/off of parameters)
      //virtual void SetThreshVal(double lwr, double upr, double off) = 0;

      virtual UnchangeableString GetName(void) = 0;
      virtual double GetTransformedVal(void) = 0;
      virtual std::string GetFixFmt(void) = 0;
      virtual double ConvertOutVal(double val) = 0;
      virtual double ConvertInVal(double val) = 0;      
      virtual const char * GetType(void) = 0;
}; /* end class ParameterABC */

/******************************************************************************
class RealParam

Represents a continuously varying parameter
******************************************************************************/
class RealParam : public ParameterABC
{
   public:      
      RealParam(void);
      RealParam(IroncladString name, double initialValue, 
                double lowerBound , double upperBound, IroncladString txIn, 
                IroncladString txOst, IroncladString txOut, IroncladString fixFmt);     
      ~RealParam(void){ DBG_PRINT("RealParam::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);

      // Declare functions for transformed values
      double GetInitialValueTransformed(void) { return m_InitialValueTransformed; };
      double GetLowerBoundTransformed(void){ return m_lowerBoundTransformed;}
      double GetUpperBoundTransformed(void){ return m_upperBoundTransformed;}
      void SetLowerBoundTransformed(double val){ m_lowerBoundTransformed = val;}
      void SetUpperBoundTransformed(double val){ m_upperBoundTransformed = val;}
      double GetEstimatedValueTransformed(void){ return m_estimatedValueTransformed;}
      double SetEstimatedValueTransformed(double estVal);

      UnchangeableString GetName(void){ return m_pName;}
      std::string GetFixFmt(void) { return m_pFixFmt; }
      double GetTransformedVal(void);
      double ConvertOutVal(double val);
      double ConvertInVal(double val);
      //threshold values (allow for implicit on/off of parameters)
      //void SetThreshVal(double lwr, double upr, double off){ m_ThreshLwr = lwr; m_ThreshUpr = upr; m_ThreshOff = off;}
      const char* GetType(void) { return "real"; };
            
   private:
      StringType m_pName;
      StringType m_pFixFmt;
      double m_InitialValueTransformed;
      double m_lowerBoundTransformed;
      double m_upperBoundTransformed;
      double m_estimatedValueTransformed;      
      //double m_ThreshLwr, m_ThreshUpr, m_ThreshOff;

      TransformTypeEnum m_TransID[NUM_STAGES];      
      void SetTransformation(TransformStageEnum which, IroncladString tx);
}; /* end class RealParam */

/******************************************************************************
class IntParam

Represents an integer parameter
******************************************************************************/
class IntParam : public ParameterABC
{
   public:      
      IntParam(void);
      IntParam(IroncladString name, int initialValue, int lowerBound, 
                   int upperBound);      
      ~IntParam(void){ DBG_PRINT("IntParam::DTOR"); Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr){sprintf(valStr, "%d", m_estimatedValueTransformed);}
      std::string GetFixFmt(void) { return m_pFixFmt; }
      double GetEstimatedValueTransformed(void){ return (double)m_estimatedValueTransformed;}
      double SetEstimatedValueTransformed(double estVal);
      double GetInitialValueTransformed(void) { return m_InitialValueTransformed; };
      double GetLowerBoundTransformed(void){ return (double)m_lowerBoundTransformed;}
      double GetUpperBoundTransformed(void){ return (double)m_upperBoundTransformed;}
      void SetLowerBoundTransformed(double val){ m_lowerBoundTransformed = (int)val;}
      void SetUpperBoundTransformed(double val){ m_upperBoundTransformed = (int)val;}
      double GetTransformedVal(void){ return (double)m_estimatedValueTransformed;}
      double ConvertOutVal(double val);
      double ConvertInVal(double val){ return val;}
      UnchangeableString GetName(void){ return m_pName;}
      void Write(FILE * pFile, int type);      
      //void SetThreshVal(double lwr, double upr, double off){ m_ThreshLwr = (int)lwr; m_ThreshUpr = (int)upr; m_ThreshOff = (int)off;}
      const char * GetType(void) {return "integer";}

   private:
      StringType m_pName;
      StringType m_pFixFmt;
      int m_InitialValueTransformed;
      int m_lowerBoundTransformed;
      int m_upperBoundTransformed;
      int m_estimatedValueTransformed;            
      //int m_ThreshLwr, m_ThreshUpr, m_ThreshOff;

      TransformTypeEnum m_TransID[NUM_STAGES];


}; /* end class IntParam */

/******************************************************************************
class ComboIntParam

Class for the combinatorial integer parameters.
******************************************************************************/
class ComboIntParam : public ParameterABC
{
   public:      
      ComboIntParam(void);
      ComboIntParam(IroncladString name, UnmoveableString configStr);
      ~ComboIntParam(void){ DBG_PRINT("ComboIntParam::DTOR"); Destroy();}
      void Destroy(void);

      void GetValAsStr(UnmoveableString valStr){sprintf(valStr, "%d", m_pCombos[m_CurIdx]);}
      void Write(FILE * pFile, int type);
      double GetLowerBoundTransformed(void){ return 0.00;}
      double GetUpperBoundTransformed(void){ return (double)(m_NumCombos - 1);}
      void SetLowerBoundTransformed(double val){ return;}
      void SetUpperBoundTransformed(double val){ return;}
      double GetEstimatedValueTransformed(void){ return (double)m_CurIdx;}
      double SetEstimatedValueTransformed(double Idx);
      UnchangeableString GetName(void){ return m_pName;}
      double GetTransformedVal(void){ return (double)(m_pCombos[m_CurIdx]);}
      double ConvertOutVal(double val){ return val;}
      double ConvertInVal(double val){ return val;}

      void SetThreshVal(double lwr, double upr, double off){ return;}
      const char * GetType(void) {return "combinatorial integer";}

   private:
      StringType m_pName;
      int m_CurIdx;
      int m_NumCombos;
      int m_InitIdx;
      int * m_pCombos;
}; /* end class ComboIntParam */

/******************************************************************************
class ComboDblParam

Class for the combinatorial real parameters.
******************************************************************************/
class ComboDblParam : public ParameterABC
{
   public:      
      ComboDblParam(void);
      ComboDblParam(IroncladString name, UnmoveableString configStr);
      ~ComboDblParam(void){ DBG_PRINT("ComboDblParam::DTOR"); Destroy();}
      void Destroy(void);

	  void GetValAsStr(UnmoveableString valStr);
      void Write(FILE * pFile, int type);
      double GetLowerBoundTransformed(void){ return 0.00;}
      double GetUpperBoundTransformed(void){ return (double)(m_NumCombos - 1);}
      void SetLowerBoundTransformed(double val){ return;}
      void SetUpperBoundTransformed(double val){ return;}
      double GetEstimatedValueTransformed(void){ return (double)m_CurIdx;}
      double SetEstimatedValueTransformed(double Idx);
      UnchangeableString GetName(void){ return m_pName;}
      double GetTransformedVal(void){ return m_pCombos[m_CurIdx];}
      double ConvertOutVal(double val){ return val;}
      double ConvertInVal(double val){ return val;}

      void SetThreshVal(double lwr, double upr, double off){ return;}
      const char * GetType(void) {return "combinatorial double";}

   private:
      StringType m_pName;
      int m_CurIdx;
      int m_NumCombos;
      int m_InitIdx;
      double * m_pCombos;      
}; /* end class ComboDblParam */

/******************************************************************************
class ComboStrParam

Class for the combinatorial string parameters.
******************************************************************************/
class ComboStrParam : public ParameterABC
{
   public:      
      ComboStrParam(void);
      ComboStrParam(IroncladString name, UnmoveableString configStr);
      ~ComboStrParam(void){ DBG_PRINT("ComboStrParam::DTOR"); Destroy();}
      void Destroy(void);

      void GetValAsStr(UnmoveableString valStr);
      void Write(FILE * pFile, int type);
      double GetLowerBoundTransformed(void){ return 0.00;}
      double GetUpperBoundTransformed(void){ return (double)(m_NumCombos - 1);}
      double GetEstimatedValueTransformed(void){ return (double)m_CurIdx;}
      void SetLowerBoundTransformed(double val){ return;}
      void SetUpperBoundTransformed(double val){ return;}
      double SetEstimatedValueTransformed(double Idx);
      UnchangeableString GetName(void){ return m_pName;}
      double GetTransformedVal(void){ return atof(m_pCombos[m_CurIdx]);}
      double ConvertOutVal(double val){ return val;}
      double ConvertInVal(double val){ return val;}

      void SetThreshVal(double lwr, double upr, double off){ return;}
      const char * GetType(void) {return "combinatorial string";}

   private:
      StringType m_pName;
      int m_CurIdx;
      int m_NumCombos;
      int m_InitIdx;
      char ** m_pCombos;
}; /* end class ComboStrParam */

/******************************************************************************
class SpecialParam

Special Ostrich parameters. These correspond to 'optimal' cost and constraint
values at any given stage of Ostrich. Can be used for linking Ostrich with
the model pre-emption capabilities of a given model.
******************************************************************************/
class SpecialParam
{
   public:      
      SpecialParam(void);
      SpecialParam(IroncladString name,  IroncladString type, 
                   IroncladString limit, IroncladString constraint, 
				   double init); 
      ~SpecialParam(void){ DBG_PRINT("SpecialParam::DTOR"); Destroy();}
      void Destroy(void);

	  void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
	  //treat special parameters as unbounded
      double GetLowerBoundTransformed(void){ return NEARLY_ZERO;}
      double GetUpperBoundTransformed(void){ return NEARLY_HUGE;}
      void SetLowerBoundTransformed(double val){ return;}
      void SetUpperBoundTransformed(double val){ return;}
      double GetEstimatedValueTransformed(void){ return m_estimatedValueTransformed;}
	  double SetEstimatedValueTransformed(double estVal){ m_estimatedValueTransformed = estVal; return 0.00;}
	  void   SetEstimatedValueTransformed(double minObj, double minCon);
	  void   SetMinObj(double minObj){ m_MinObj = minObj;}
      UnchangeableString GetName(void){ return m_pName;}
      double GetTransformedVal(void){ return m_estimatedValueTransformed;}
      double ConvertOutVal(double val){ return val;}
      double ConvertInVal(double val){ return val;}
	  void   Enable(void){ m_bSet = true;}
      //threshold values (allow for implicit on/off of parameters)
      void SetThreshVal(double lwr, double upr, double off){ return;}
      const char * GetType(void) {return "special";}
	  ConstraintABC * GetConstraint(void);
            
   private:
      StringType m_pName;
	  StringType m_pType; 
      StringType m_pLimit;
	  StringType m_pConstraint;
	  double m_MinObj;
      double m_estimatedValueTransformed;
	  bool m_bSet;
}; /* end class SpecialParam */


#endif /* PARAMETER_ABC_H */
