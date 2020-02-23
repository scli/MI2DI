#ifndef _MIC_PARA_MGR_
#define _MIC_PARA_MGR_

enum DataOpt {_THE_21_, _EXPRESSION_, _OTHERS_};
enum InputType { _PROTEIN_, _DNA_};
enum InvOpt  {_LASSO_, _REMOVE_ONE_, _REGULARIZATION_};
enum MatrixOpt { _SANDER_, _LOCAL_WEIGHT_};
enum PartitionOpt { _UNIFORM_PARTITION_, _TO_BE_DECIDED_};
enum InvCorrOpt   {_INV_COV_MAT_FIRST_, _INV_CORR_MAT_FIRST_};

class MICParaMgr
{
   public:
     MICParaMgr();
     ~MICParaMgr();
     static bool         _ALIGN_SKIP_END_;
     static int          _ALIGN_SKIP_END_COUNTER_;
     static double       _ALIGN_FILTER_END_THRES_;
     static double       _ALIGN_FILTER_COL_THRES_;
     static DataOpt      _DATA_OPT_;
     
     static InputType    _INPUT_TYPE_;
     static InvOpt       _INVERSE_OPT_;
     static MatrixOpt    _MATRIX_MODE_;
     static PartitionOpt _DATA_PARTITION_OPT_;
     static InvCorrOpt   _INV_CORR_OPT_;
     static int          _NUM_PARTITION_CLUSTERS_;
     static double       _PSEUDO_COUNT_WEIGHT_;


};


#endif
