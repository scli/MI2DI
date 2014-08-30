#include "MICParaMgr.h"
bool         MICParaMgr::_ALIGN_SKIP_END_         = false;
double       MICParaMgr::_ALIGN_FILTER_END_THRES_ = 0.2;
int          MICParaMgr::_ALIGN_SKIP_END_COUNTER_ = 3;
double       MICParaMgr::_ALIGN_FILTER_COL_THRES_ =0.4;
DataOpt      MICParaMgr::_DATA_OPT_               =_THE_21_;
InputType    MICParaMgr::_INPUT_TYPE_             =_PROTEIN_; 
//InvOpt    MICParaMgr::_INVERSE_OPT_            =_REGULARIZATION_;
InvOpt       MICParaMgr::_INVERSE_OPT_            =_REMOVE_ONE_;
MatrixOpt    MICParaMgr::_MATRIX_MODE_            =_SANDER_;
//MatrixOpt MICParaMgr::_MATRIX_MODE_            =_LOCAL_WEIGHT_;
PartitionOpt MICParaMgr::_DATA_PARTITION_OPT_     =_UNIFORM_PARTITION_;
InvCorrOpt   MICParaMgr::_INV_CORR_OPT_           = _INV_COV_MAT_FIRST_;
int          MICParaMgr::_NUM_PARTITION_CLUSTERS_ =3;
double       MICParaMgr::_PSEUDO_COUNT_WEIGHT_    =0.5;
MICParaMgr::MICParaMgr()
{
   



}
