#ifndef _AGGREGATION_H_
#define _AGGREGATION_H_

#include "MICData.h"
#include <vector>
using namespace std;
struct ScoreItem
{
   int i1;
   int i2;
   float score;
};
class Aggregation
{
    protected:    
        MICData* mData;
        double** mMIMatrix;

        void getMatrixW(int pos1, int pos2, double** w);
        
        void vectorTimesMatrix(double* vect, double** mat, int dim1, int dim2, double * result);
        void vectorTimesTranspose(double* vect, double** mat, int dim1, int dim2, double * result);
        void normalize(double* vect, int size);
        void normalize(double** vect, int dim1, int dim2);
        void innerDivide(double* v1, double* v2, double* res, int size);
        double maximum(double* v1, double* v2, int size);

        void topK(ScoreItem * items, int start, int end, int k, double pivot);

    public:
        Aggregation(MICData* data);
        ~Aggregation();
        void topMIPairs(vector<int>* offset);


};
#endif
