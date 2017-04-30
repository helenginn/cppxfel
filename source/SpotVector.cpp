//
//  SpotVector.cpp
//  cppxfel
//
//  Created by Helen Ginn on 15/10/2015.
//  Copyright (c) 2015 Division of Structural Biology Oxford. All rights reserved.
//

#include "SpotVector.h"
#include "Matrix.h"
#include "misc.h"
#include "FileParser.h"
#include "Detector.h"
#include "PNGFile.h"

double SpotVector::trustComparedToStandardVector(SpotVectorPtr standardVector)
{
    double standardDistance = standardVector->distance();
    double diff = fabs(standardDistance - distance());
    
    return 1 / diff;
}

SpotVector::SpotVector(vec transformedHKL, vec normalHKL)
{
    firstSpot = SpotPtr();
    secondSpot = SpotPtr();
    approxResolution = 0;
    minDistanceTolerance = 0;
	_isIntraPanelVector = -1;

    update = false;
    hkl = normalHKL;
    spotDiff = copy_vector(transformedHKL);
    calculateUnitVector();
}

SpotVector::SpotVector(SpotPtr first, SpotPtr second)
{
    firstSpot = first;
    secondSpot = second;
    update = false;
    approxResolution = 0;
    minDistanceTolerance = 0;
    minAngleTolerance = 0;
	_isIntraPanelVector = -1;

    if (!first || !second)
        return;
    
    calculateDistance();
    firstDistance = cachedDistance;
}

void SpotVector::calculateUnitVector()
{
    cachedDistance = length_of_vector(spotDiff);
    unitSpotDiff = copy_vector(spotDiff);
    scale_vector_to_distance(&unitSpotDiff, 1.);
}

void SpotVector::quickDistance()
{
	if (firstSpot && secondSpot)
	{
		vec firstVector = firstSpot->storedVector();
		vec secondVector = secondSpot->storedVector();

		spotDiff = copy_vector(secondVector);
		take_vector_away_from_vector(firstVector, &spotDiff);
	}

	calculateUnitVector();
}

void SpotVector::calculateDistance()
{
    if (firstSpot && secondSpot)
    {
        vec firstVector = firstSpot->estimatedVector();
        vec secondVector = secondSpot->estimatedVector();
        
        spotDiff = copy_vector(secondVector);
        take_vector_away_from_vector(firstVector, &spotDiff);
    }
    
    calculateUnitVector();
}

double SpotVector::distance()
{
    if (update || cachedDistance != cachedDistance)
    {
        calculateDistance();
        update = false;
    }
    
    return cachedDistance;
}

// in radians
double SpotVector::angleWithVector(SpotVectorPtr spotVector2)
{
    return angleBetweenVectors(spotVector2->getUnitVector(), unitSpotDiff, true);
}

double SpotVector::cosineWithVector(SpotVectorPtr spotVector2)
{
    return cosineBetweenUnitVectors(spotVector2->getUnitVector(), unitSpotDiff);
}

double SpotVector::cosineWithVertical()
{
    return cosineBetweenUnitVectors(getUnitVector(), new_vector(0, 1, 0));
}

bool SpotVector::isCloseToSpotVector(SpotVectorPtr spotVector2, double maxDistance)
{
    vec spotDiff2 = spotVector2->getVector();
    
    if (fabs(spotDiff2.h - spotDiff.h) > maxDistance || fabs(spotDiff2.k - spotDiff.k) > maxDistance
        || fabs(spotDiff2.l - spotDiff.l) > maxDistance)
        return false;
    
    return true;
}

double SpotVector::similarityToSpotVector(SpotVectorPtr spotVector2)
{
    vec displace = copy_vector(spotDiff);
    take_vector_away_from_vector(spotVector2->getVector(), &displace);
    
    return length_of_vector(displace);
}

SpotVectorPtr SpotVector::copy()
{
    SpotVectorPtr newPtr = SpotVectorPtr(new SpotVector(firstSpot, secondSpot));
    newPtr->hkl = copy_vector(hkl);
    newPtr->spotDiff = copy_vector(spotDiff);
    newPtr->update = update;
    newPtr->sameLengthStandardVectors = sameLengthStandardVectors;
    
    return newPtr;
}

SpotVectorPtr SpotVector::vectorRotatedByMatrix(MatrixPtr mat)
{
    SpotVectorPtr newVec = copy();
    
    mat->multiplyVector(&newVec->hkl);
    mat->multiplyVector(&newVec->spotDiff);
    
    return newVec;
}

bool SpotVector::hasCommonSpotWithVector(SpotVectorPtr spotVector2)
{
    if (firstSpot == spotVector2->firstSpot || firstSpot == spotVector2->secondSpot)
    {
        return true;
    }
    
    if (secondSpot == spotVector2->firstSpot || secondSpot == spotVector2->secondSpot)
    {
        return true;
    }
    
    return false;
}

std::string SpotVector::description()
{
    return "(" + f_to_str(spotDiff.h) + ", " + f_to_str(spotDiff.k) + ", " + f_to_str(spotDiff.l) + ")";
}

void SpotVector::addSimilarLengthStandardVectors(std::vector<SpotVectorPtr> standardVectors, double tolerance)
{
    sameLengthStandardVectors.clear();
    tolerance = this->getMinDistanceTolerance();
    
    for (int i = 0; i < standardVectors.size(); i++)
    {
        double trust = trustComparedToStandardVector(standardVectors[i]);
        
        if (trust > tolerance)
        {
            sameLengthStandardVectors.push_back(standardVectors[i]);
        }
    }
    
    std::ostringstream logged;
    logged << "Added " << sameLengthStandardVectors.size() << " similar standard lengths to vector (min tolerance " << getMinDistanceTolerance() << ")." << std::endl;
    Logger::mainLogger->addStream(&logged, LogLevelDebug);
}

SpotVectorPtr SpotVector::differenceFromVector(SpotVectorPtr spotVec)
{
    vec spotDiffDiff = copy_vector(spotDiff);
    
    take_vector_away_from_vector(spotVec->spotDiff, &spotDiffDiff);
    
    SpotVectorPtr newVec = SpotVectorPtr(new SpotVector(spotDiffDiff, new_vector(0, 0, 0)));
    
    return newVec;
}

SpotVectorPtr SpotVector::vectorBetweenSpotsFromArray(std::vector<SpotVectorPtr> vectors, SpotPtr spot1, SpotPtr spot2)
{
    for (int i = 0; i < vectors.size(); i++)
    {
        SpotPtr firstSpot = vectors[i]->getFirstSpot();
        
        if (firstSpot == spot1 || firstSpot == spot2)
        {
            SpotPtr secondSpot = vectors[i]->getSecondSpot();
            
            if ((secondSpot == spot1 || secondSpot == spot2) && secondSpot != firstSpot)
            {
                return vectors[i];
            }
        }
    }
    
    return SpotVectorPtr();
}

double SpotVector::getResolution()
{
    if (approxResolution != 0)
        return approxResolution;
    
    double resol1 = firstSpot->resolution();
    double resol2 = secondSpot->resolution();
    
    double minResolution = std::max(resol1, resol2);
    
    approxResolution = minResolution;
    
    return approxResolution;
}

double SpotVector::getMinDistanceTolerance()
{
    if (minDistanceTolerance != 0)
        return minDistanceTolerance;
    
    double resolution = getResolution();
    
    double rlpSize = FileParser::getKey("INDEXING_RLP_SIZE", 0.001);
    
    //double mosaicity = FileParser::getKey("INITIAL_MOSAICITY", 0.0);
    double bandwidth = FileParser::getKey("INITIAL_BANDWIDTH", 0.0013);
    double wavelength = FileParser::getKey("INTEGRATION_WAVELENGTH", 0.);
    
    double minWavelength = wavelength * (1 - bandwidth * 2);
    double maxWavelength = wavelength * (1 + bandwidth * 2);
    
    // we set k = 0
    
    double minRadius = 1 / minWavelength;
    double minL = - pow(resolution, 2) / (2 * minRadius);
    double minH = sqrt(pow(resolution, 2) - pow(minL, 2));
    
    double maxRadius = 1 / maxWavelength;
    double maxL = - pow(resolution, 2) / (2 * maxRadius);
    double maxH = sqrt(pow(resolution, 2) - pow(maxL, 2));
    
    vec minVec = new_vector(minH, 0, minL);
    vec maxVec = new_vector(maxH, 0, maxL);
    
    take_vector_away_from_vector(minVec, &maxVec);
    
    minDistanceTolerance = length_of_vector(maxVec);
    minDistanceTolerance += rlpSize;
    
    minDistanceTolerance = 1 / minDistanceTolerance;
    
    return minDistanceTolerance;
}

bool SpotVector::isIntraPanelVector()
{
	if (_isIntraPanelVector >= 0)
	{
		return _isIntraPanelVector;
	}

    DetectorPtr onePanel = firstSpot->getDetector();
    DetectorPtr twoPanel = secondSpot->getDetector();
    
    if (!onePanel || !twoPanel)
	{
		_isIntraPanelVector = 0;
        return _isIntraPanelVector;
	}

    if ((onePanel != twoPanel) && onePanel->getParent() == twoPanel->getParent())
    {
        if (!onePanel->getParent()->isRefinable(GeometryScoreTypeInterpanel))
        {
			_isIntraPanelVector = 1;
            return _isIntraPanelVector;
        }
    }

	_isIntraPanelVector = (onePanel == twoPanel);
    return _isIntraPanelVector;
}

bool SpotVector::spansChildrenOfDetector(DetectorPtr parent)
{
    DetectorPtr onePanel = firstSpot->getDetector();
    DetectorPtr twoPanel = secondSpot->getDetector();
    
    if (!onePanel || !twoPanel)
    {
        return false;
    }
    
    if (onePanel == twoPanel)
    {
        return false;
    }
    
    if (parent->isAncestorOf(onePanel) && parent->isAncestorOf(twoPanel))
    {

        for (int i = 0; i < parent->childrenCount(); i++)
        {
            if (parent->getChild(i)->isAncestorOf(onePanel) &&
                parent->getChild(i)->isAncestorOf(twoPanel))
            {
                return false;
            }
        }
        
        return true;
    }
    
    return false;
}

bool SpotVector::isOnlyFromDetector(DetectorPtr detector)
{
    DetectorPtr onePanel = firstSpot->getDetector();
    DetectorPtr twoPanel = secondSpot->getDetector();
    
    return (detector->isAncestorOf(onePanel) && detector->isAncestorOf(twoPanel));
}

double SpotVector::getMinAngleTolerance()
{
    if (minAngleTolerance != 0)
    {
        return minAngleTolerance;
    }
    
    double rlpSize = FileParser::getKey("INDEXING_RLP_SIZE", 0.001);
    double distTolerance = getMinDistanceTolerance();
    double realDistance = distance();
    
    double reciprocalFatnessRadius = (rlpSize * 2 + 1 / distTolerance) / 2;
    
    double tanAngle = tan(reciprocalFatnessRadius / realDistance);
    double angle = atan(tanAngle) * 2;
    
    minAngleTolerance = angle;
    
    return minAngleTolerance;
}

bool SpotVector::usesBeamCentre()
{
    return (firstSpot->isBeamCentre() || secondSpot->isBeamCentre());
}

bool SpotVector::originalDistanceLessThan(double threshold)
{
    return (firstDistance < threshold);
}

void SpotVector::drawOnImage(PNGFilePtr file)
{
    Coord spot1 = getFirstSpot()->getXY();
    Coord spot2 = getSecondSpot()->getXY();
    
    file->drawLine(spot1.first, spot1.second, spot2.first, spot2.second,
                   0.0, 0, 0, 100);
}
