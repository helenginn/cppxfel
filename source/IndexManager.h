//
//  IndexManager.h
//  cppxfel
//
//  Created by Helen Ginn on 14/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef __cppxfel__IndexManager__
#define __cppxfel__IndexManager__

#include "cmtzlib.h"
#include "csymlib.h"
#include <stdio.h>
#include <vector>
#include <map>
#include "Image.h"
#include "Vector.h"
#include "parameters.h"
#include "Logger.h"

class IndexManager
{
protected:
    std::vector<Image *> images;
    std::vector<double> unitCell;
    MatrixPtr unitCellOnly;
    MatrixPtr unitCellMatrix;
    MatrixPtr unitCellMatrixInverse;
    CSym::CCP4SPG *spaceGroup;
    int spaceGroupNum;
    std::vector<MtzPtr> mtzs;
    double minimumTrustDistance;
    double minimumTrustAngle;
    double solutionAngleSpread;
    
    double maxDistance;
    std::ostringstream logged;
    std::vector<VectorDistance> vectorDistances;
public:
    Image *getImage(int i)
    {
        return images[i];
    }
    
    std::vector<MtzPtr> getMtzs()
    {
        return mtzs;
    }
    
    static void indexThread(IndexManager *indexer, std::vector<MtzPtr> *mtzSubset, int offset);
    void index();
    IndexManager(std::vector<Image *>images);
    
    void sendLog(LogLevel priority = LogLevelNormal);
};

#endif /* defined(__cppxfel__IndexManager__) */
