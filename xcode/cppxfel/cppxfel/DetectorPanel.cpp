//
//  DetectorPanel.cpp
//  cppxfel
//
//  Created by Helen Mary Elizabeth Duyvesteyn on 05/11/2015.
//  Copyright © 2015 Structural Biology. All rights reserved.
//

#include "DetectorPanel.hpp"
#include "Panel.h"
#include "Miller.h"
#include "Logger.h"

#import "misc.h" //Conversion of degress/radians.
#import "FileParser.h"

//MM_PER_PIXEL

//Declare memory pointer.
DetectorPanelPtr DetectorPanel::masterPanel = DetectorPanelPtr();

//Constructor for detector panel class.
//Assign height and width in header.

DetectorPanel::DetectorPanel (double new_height, double new_width, double new_rot_x, double new_rot_y, double new_rot_z, double new_centre_x, double new_centre_y, double new_centre_z)
{
    height = new_height;
    width = new_width;
    rot_x = deg2rad(new_rot_x);
    rot_y = deg2rad(new_rot_y);
    rot_z = deg2rad(new_rot_z);
    //Corresponds to beam centre.
    centre_x = new_centre_x;
    centre_y =new_centre_y;
    centre_z = new_centre_z;
    
    
    //Parse in degrees. Conversion in misc.h.
    //Set static variables to initial values.

}

//Getters and setters. get_masterpanel. Get one detector's height and width.
//HEADER FILE. PRIVATE variables.

//THE MASTER.
DetectorPanelPtr DetectorPanel::getMasterPanel()
{
    if (!masterPanel) //False ==empty
    {
        setupMasterPanel();
    }
    return masterPanel;
}
//Panel up from wherever.



void DetectorPanel::setSuperPanel(boost::weak_ptr<DetectorPanel> toMotherPtr)
{
    superPanel = toMotherPtr;
}
//Weak pointer constructor, such that little panels can be told what the big panel is called.
/////////////////////////////////


//Within subpanels vector, create a newPanel pointer/array which calls on superpanel and constructs a newPanel for each subpanel. Append to subpanel array using addSubPanel function.

DetectorPanelPtr DetectorPanel::getSubPanel(int i)
{
    return subPanels[i];
}

void DetectorPanel::createSubPanel(double topLeft_x, double topLeft_y, double bottomRight_x, double bottomRight_y, DetectorPanelPtr superPanel) //Vector should contain newPanels.
{
    double mmPerPixel = FileParser::getKey("MM_PER_PIXEL", 0.11); //Default for cspad.
    double centre_x = mmPerPixel*((bottomRight_x-topLeft_x)/2);
    double centre_y = mmPerPixel*(bottomRight_y-topLeft_y)/2;
    
    double uppervec_y = (mmPerPixel * topLeft_y) - height/2;
    double uppervec_x = (mmPerPixel*topLeft_x) - width/2;
    double actual_centre_y = uppervec_y + centre_y;
    double actual_centre_x = uppervec_x + centre_x;
    double actual_centre_z = 0;
    
    double new_width = centre_x * 2;
    double new_height = centre_y * 2;
    double rot_x = 0;
    double rot_y = 0;
    double rot_z = 0;
    
    DetectorPanelPtr panel = DetectorPanelPtr (new DetectorPanel(new_height, new_width, rot_x, rot_y, rot_z, actual_centre_x, actual_centre_y, actual_centre_z));
    
    //newPanel is a new subPanel. Add to subPanel vector.
    panel->setSuperPanel(selfPtr); //Return superPanel Ptr to original.
    addSubPanel(panel); //Add subpanel to subpanels vector.
}
void DetectorPanel::addSubPanel(DetectorPanelPtr newPanel)//Adds newPanel to subPanel matrix.
{
    subPanels.push_back(newPanel);
    newPanel->setSuperPanel(selfPtr);
}
//actualcentre z = 0

/////////////////////////////////////
void DetectorPanel::setupMasterPanel()
{
    double mmPerPixel = FileParser::getKey("MM_PER_PIXEL", 0.11); //Default for cspad.
    double detectorDistance = FileParser::getKey("DETECTOR_DISTANCE", 0.0);
    std::vector<int> imageSize = FileParser::getKey("DETECTOR_SIZE", std::vector<int>());
    double delta_x = 0.0;
    double delta_y = 0.0;
    double delta_z = 0.0;
    
    if (imageSize.size() != 2)
    {
        std::cout << "Setting to default image dimensions for CSPAD." << std::endl;
        imageSize.clear();
        imageSize.push_back(1765);
        imageSize.push_back(1765);
    }
    
    std::vector<double> beamCentre = FileParser::getKey("BEAM_CENTRE", std::vector<double>());
    if (beamCentre.size() != 2)
    {
        std::cout << "Setting to default beam centre." << std::endl;
        beamCentre.clear();
        beamCentre.push_back (imageSize[0] / 2);
        beamCentre.push_back (imageSize[1] / 2);
    }
    
    double actualCentre_x = (imageSize[0] / 2);
    double actualCentre_y = (imageSize[1] / 2);
    
    double diff_Centre_x = beamCentre[0] - actualCentre_x;
    double diff_Centre_y = beamCentre[1] - actualCentre_y;
    double rot_y = FileParser::getKey("DETECTOR_ROTATION_Y", 0.0);
    
    DetectorPanel::convertPixelstoOffset(rot_y, diff_Centre_x, diff_Centre_y, mmPerPixel, &delta_x, &delta_y, &delta_z);
    delta_z += detectorDistance;
    
    double new_height = imageSize[1] * mmPerPixel;
    double new_width = imageSize[0] * mmPerPixel;
    
    masterPanel = DetectorPanelPtr(new DetectorPanel (new_height, new_width, rot_y, 0.0, 0.0, delta_x, delta_y, delta_z));
    
}

void DetectorPanel::convertPixelstoOffset(double rot_y, double diff_x, double diff_y, double mmPerPixel, double *delta_x, double *delta_y, double *delta_z)
{
    rot_y = deg2rad(rot_y);
    *delta_z = (sin (rot_y) * diff_x);
    *delta_x = (cos (rot_y) * diff_x);
    *delta_y = diff_y; //CHANGE.
    
    *delta_z *= mmPerPixel;
    *delta_x *= mmPerPixel;
    *delta_y *= mmPerPixel;
    
}


    //Convert panel pixels into mm THEN sinx, cosy.
    
    
    //Template for pixel number.
//Turn user input into panels.

//miller2xyz (miller, coord)


//Description of properties of panel...std::cout.

void DetectorPanel::descriptionPanel()
{
    std::ostringstream logged;
    logged << "Setting panel rotations to " << rot_x << ", " << rot_y << ", " << rot_z << std::endl;
    
    logged << "Setting panel origin to " << centre_x << ", " << centre_y << ", " << centre_z << std::endl;
    Logger::mainLogger->addStream(&logged);
}


//Subpanel inherits origin-->relative origin.


//Temp. Use with detector panels??? Use vec origin.



//cf Miller::miller(ARGUMENTS) in miller.cpp.

//Implement detector panel pointer. *
//std::vector detector panel pointers.
//Every panel. Define origin in centre. Double centre x/y
//Three angles. Rotx/Roty/Rotz.
//Parse in parameters.

//MM_PER_PIXEL variable associated with detector panel. Get key from file parser.-INPUT. Should not have to enter into constructor.
//Recursive functions to set origin. Origin asks first superpanel what origin is. Add x,y,z onto a vector.

//If x% panel is negative background then discount panel.


