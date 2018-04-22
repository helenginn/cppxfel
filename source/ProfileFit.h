//
//  ProfileFit.h
//  cppxfel
//
//  Created by Helen Mary Elizabeth Duyvesteyn on 04/03/2018.
//  Copyright © 2018 Structural Biology. All rights reserved.
//

#ifndef PROFILEFIT_H_
#define PROFILEFIT_H_

#include <stdio.h>
#include <string>
#include <map>
#include "parameters.h"
#include "Shoebox.h"
#include "Miller.h"


class ProfileFit : public Shoebox
{
private:
    //int shoeBoxF;
    //int shoeBoxB;
    //int shoeBoxN;
    
protected:
    std::string filename;
    
public:
   // ~ProfileFit();
    //Profile fits for image???
    void calculateImageProfile(ImagePtr image);
    void addShoeboxes(ShoeboxPtr shoeTmp);
    void subtractfromShoebox(double toSubtract);
    ProfileFit() : Shoebox(MillerPtr()){};
    void correctSpotCentring(double *xCoord, double *yCoord, ImagePtr img);
    ShoeboxPtr intensitySearchGrid(int startX, int startY, ImagePtr im, std::string fName, int spotNum, double *bg);
    double estimateFWHM(double stdevI = 0);
    void makeShoebox(ShoeboxPtr); //Not used.
};




#endif /* PROFILEFIT_H_ */
