//
//  DetectorPanel.hpp
//  cppxfel
//
//  Created by Helen Mary Elizabeth Duyvesteyn on 05/11/2015.
//  Copyright © 2015 Structural Biology. All rights reserved.
//

#ifndef DetectorPanel_hpp
#define DetectorPanel_hpp
#import <vector>
#import "parameters.h"
#import <boost/weak_ptr.hpp>
#include <stdio.h>
#include <math.h>
#include<cmath>
#import "FileParser.h"


//POINTERS std::shared_pointer
//DetectorPanel Pointer.
//Parameters.

class DetectorPanel
{
    //origin (double centre_x, centre_y);
private:
    //Master Panel
    double height;
    double width;
    
    double rot_x;
    double rot_y;
    double rot_z;
    
    double centre_x;
    double centre_y;
    double centre_z;
    
    double pixel_beamCentre_x;//USED?
    double pixel_beamCentre_y;//USED?
    double mmPerPixel;
    
    static DetectorPanelPtr masterPanel; //Keeps reference, boost share pointer.
    
    std::vector<DetectorPanelPtr> subPanels; //TREEE. Array with pointers rather than obj.
    
    boost::weak_ptr<DetectorPanel> superPanel;
    
    boost::weak_ptr<DetectorPanel> selfPtr;
    
    static void setupMasterPanel();
    
    
    double rotX;
    
    
    
public:
    //getters. In cpp, uses function as actual function. If declare in header, compiled into whatever calling. Improves efficiency, but increases file.
    double getHeight()
    {
        return height;
    }
    double getWidth()
    {
        return width;
    }
    
    
    static DetectorPanelPtr getMasterPanel();
    //DetectorPanelPtr getSuperPanel();
    void setSuperPanel(boost::weak_ptr<DetectorPanel> toMotherPtr);
    
    static void convertPixelstoOffset(double rot_x, double diff_x, double diff_y, double mmPerPixel, double *delta_x, double *delta_y, double *delta_z);
    
    //Constructor, ergo no type.
    DetectorPanel (double new_height, double new_width, double new_rot_x, double new_rot_y, double new_rot_z, double new_centre_x, double new_centre_y, double new_centre_z);
    DetectorPanelPtr getSubPanel(int i);
    
    //Declare constructor in header.
    
    DetectorPanel *parentDetectorPanel;
    //Getter static for MasterPanel
    
    void descriptionPanel();
    
    void createSubPanel(double topLeft_x, double topLeft_y, double bottomRight_x, double bottomRight_y, DetectorPanelPtr superPanel);
    
    ////////////////////////////////
    ///SUBPANELS ETC//////////////
    //////////////////////////////
    
    
    //To do: Create getters for x,y,z, Rotx.
    //Extract VALUE *& pointed to by weak ptr superPanel and ADD to input values.
    //Return above.
    
    //getrot_x (centre_x,y,z) and returns its own rotation...sincostan.
    //so getRotX() should call the superpanel's getRotX() and add on its result
    
    //
    double bottomRight_x;
    double topLeft_x;
    
    
    double pythagoras (double x1, double x2, double y1, double y2) //ADD to cMisc?
    {
        double hypotenuse = sqrt(pow(x2-x1,2) + pow(y2-y1,2));
        return hypotenuse;
    }
    
    
    void addSubPanel(DetectorPanelPtr newPanel); //subPanels in private is array with all pointers.
    
    
    DetectorPanelPtr getSuperPanel() //Provides weak pointer with temporary access to DetectorPanelPtr.
    {
        return superPanel.lock();
    }
    
    
    
    
    double getRotY()
    {
        if (!getSuperPanel()) //If no superPanel value, just origin of 0.
        {
            return 0;
        }
        
        double new_rot_y = getSuperPanel()->getRotY();
        new_rot_y += rot_y;
        return new_rot_y;
        
    }
    
    //deg2rot?
    
    
    
    
    //x,y user input parameters.
    
    
    double getCentreX()
    {
        if (!getSuperPanel()) //If no superPanel value, just origin of 0.
        {
            return 0;
        }
        double new_centre_x = getSuperPanel()->getCentreX();
        //Find address of new_centre_x from detector panel. CREATE CENTREX().
        
        double rot_y_correction = getSuperPanel()->getRotY();
        double hypotenuse = pythagoras(centre_x, 0, centre_y, 0);
        new_centre_x += hypotenuse * acos(rot_y_correction);
        //Consider other rotations.
        new_centre_x += centre_x;
        
        //Use value of new_centre_x to find offset for centre_x.
        //Initialise sub_x, or similar.
        return new_centre_x; //Modified centrex. If rotation in masterpanel, then alter coords.
    }
    double getCentreY()
    {
        if (!getSuperPanel()) //If no superPanel value, just origin of 0.
        {
            return 0;
        }
        double new_centre_y = getSuperPanel()->getCentreY();
        //double new_rot_x = getSuperPanel()->getRotX();
        new_centre_y += centre_y;
        return centre_y;
    }
    double getCentreZ()
    {
        
        if (!getSuperPanel()) //If no superPanel value, just origin of 0.
        {
            return 0;
        }
        double new_centre_z = getSuperPanel()->getCentreZ();
        double rot_y_correction = getSuperPanel()->getRotY();
        double hypotenuse = pythagoras(centre_x, 0, centre_y, 0);
        new_centre_z += hypotenuse * asin(rot_y_correction);
        //Consider other rotations.
        new_centre_z += centre_z;
        return centre_z;
    }
    
    ///////////////
    
    //When declare detector panel, get pointer to detector panel. Put into object. Tells to take a weak reference.
    
};

//Getters and setters.

#endif