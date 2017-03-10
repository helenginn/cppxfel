/*
 * Vector.c
 *
 *  Created on: 27 Aug 2014
 *      Author: helenginn
 */

#include "Vector.h"
#include <cmath>
#include <cfloat>
#include <iostream>
#include <tuple>
#include "Matrix.h"
#include "Logger.h"
#include <fstream>



vec reverseVector(vec vec1)
{
    vec vec2 = copy_vector(vec1);
    vec2.h = -vec2.h;
    vec2.k = -vec2.k;
    vec2.l = -vec2.l;
    
    return vec2;
}

bool vectors_are_equal(vec vec1, vec vec2)
{
    return (vec1.h == vec2.h && vec1.k == vec2.k && vec1.l == vec2.l);
}

double distance_between_vectors(vec vec1, vec vec2)
{
    take_vector_away_from_vector(vec1, &vec2);
    
    return length_of_vector(vec2);
}

bool within_vicinity(vec vec1, vec vec2, double maxD)
{
    if (fabs(vec1.h - vec2.h) > maxD)
        return false;
    
    if (fabs(vec1.k - vec2.k) > maxD)
        return false;
    
    if (fabs(vec1.l - vec2.l) > maxD)
        return false;
    
    take_vector_away_from_vector(vec1, &vec2);
    
    if (length_of_vector(vec2) > maxD)
        return false;
    
    return true;
}

vec vector_between_vectors(vec vec1, vec vec2)
{
	vec vec;
	vec.h = vec2.h - vec1.h;
	vec.k = vec2.k - vec1.k;
	vec.l = vec2.l - vec1.l;

	return vec;
}

void take_vector_away_from_vector(vec vec1, vec *vec2)
{
	(*vec2).h -= (vec1).h;
	(*vec2).k -= (vec1).k;
	(*vec2).l -= (vec1).l;
}

void add_vector_to_vector(vec *vec1, vec vec2)
{
	(*vec1).h += vec2.h;
	(*vec1).k += vec2.k;
	(*vec1).l += vec2.l;
}

MatrixPtr rotation_between_vectors(vec vec1, vec vec2)
{
    MatrixPtr matrix = MatrixPtr(new Matrix());
   // MatrixPtr matrix;
    
    // Find closest angle between vectors (straightest sweep between the two)
    double cosine = cosineBetweenVectors(vec1, vec2);
    
    vec crossVector = cross_product_for_vectors(vec1, vec2);
    scale_vector_to_distance(&crossVector, 1);
    
    double angle = acos(cosine);
    
    matrix->rotateRoundUnitVector(crossVector, angle);
    
    return matrix;
}

MatrixPtr closest_rotmat_analytical(vec vec1, vec vec2,
                                    vec axis, double *resultantAngle, bool addPi)
{
    scale_vector_to_distance(&vec1, 1);
    scale_vector_to_distance(&vec2, 1);
    scale_vector_to_distance(&axis, 1);
    
    double a = vec2.h;
    double b = vec2.k;
    double c = vec2.l;
    
    double p = vec1.h;
    double q = vec1.k;
    double r = vec1.l;
    
    double x = axis.h;
    double y = axis.k;
    double z = axis.l;
    
    double A = a*(p*x*x - p + x*y*q + x*z*r) +
    b*(p*x*y + q*y*y - q + r*y*z) +
    c*(p*x*z + q*y*z + r*z*z - r);
    
    double B = a*(y*r - z*q) + b*(p*z - r*x) + c*(q*x - p*y);
    
    double tanTheta = - B / A;
    double theta = atan(tanTheta);
    
   // work out why this is always correct...? NOPE IT ISN'T YOU IDIOT
    // ARRRRRGH
    
//    *resultantAngle = bestAngle;
    double cc = cos(theta);
    double C = 1 - cc;
    double s = sin(theta);
    double occ = -cc;
    double oC = 1 - occ;
    double os = -s;
    
    double pPrime = (x*x*C+cc)*p + (x*y*C-z*s)*q + (x*z*C+y*s)*r;
    double qPrime = (y*x*C+z*s)*p + (y*y*C+cc)*q + (y*z*C-x*s)*r;
    double rPrime = (z*x*C-y*s)*p + (z*y*C+x*s)*q + (z*z*C+cc)*r;
    
    double pDoublePrime = (x*x*oC+occ)*p + (x*y*oC-z*os)*q + (x*z*oC+y*os)*r;
    double qDoublePrime = (y*x*oC+z*os)*p + (y*y*oC+occ)*q + (y*z*oC-x*os)*r;
    double rDoublePrime = (z*x*oC-y*os)*p + (z*y*oC+x*os)*q + (z*z*oC+occ)*r;
    
    double cosAlpha = pPrime * a + qPrime * b + rPrime * c;
    double cosAlphaOther = pDoublePrime * a + qDoublePrime * b + rDoublePrime * c;
    
    addPi = (cosAlphaOther > cosAlpha);
    double bestAngle = theta + addPi * M_PI;
    
    MatrixPtr mat = MatrixPtr(new Matrix());
    mat->rotateRoundUnitVector(axis, bestAngle);
    
    return mat;
}

vec cross_product_for_vectors(vec vec1, vec vec2)
{
	double new_h = vec1.k * vec2.l - vec1.l * vec2.k;
	double new_k = vec1.l * vec2.h - vec1.h * vec2.l;
	double new_l = vec1.h * vec2.k - vec1.k * vec2.h;

	return new_vector(new_h, new_k, new_l);
}

double dot_product_for_vectors(vec vec1, vec vec2)
{
    return vec1.h * vec2.h + vec1.k * vec2.k + vec1.l * vec2.l;
}

vec copy_vector(vec old_vec)
{
	vec vec;
    memcpy(&vec, &old_vec, sizeof(vec.h) * 3);

	return vec;
}

double cosineBetweenUnitVectors(vec vec1, vec vec2)
{
    double dotProduct = vec1.h * vec2.h + vec1.k * vec2.k + vec1.l * vec2.l;
    
    return dotProduct;
}

double cosineBetweenVectors(vec vec1, vec vec2)
{
    double dotProduct = vec1.h * vec2.h + vec1.k * vec2.k + vec1.l * vec2.l;
    double vec1_length = length_of_vector(vec1);
    double vec2_length = length_of_vector(vec2);
    
    double cosTheta = dotProduct / (vec1_length * vec2_length);
    
    return cosTheta;
}

double angleBetweenVectors(vec vec1, vec vec2, int isUnit)
{
    double cosTheta = 0;
    
    if (!isUnit)
    {
        cosTheta = cosineBetweenVectors(vec1, vec2);
    }
    else
    {
        cosTheta = cosineBetweenUnitVectors(vec1, vec2);
    }

    if (cosTheta > 1)
        cosTheta = 1;
    
    if (cosTheta < -1)
        cosTheta = -1;
    
	double angle = acos(cosTheta);

    if (angle != angle && (cosTheta < 1.0001))
        angle = 0;
    
	return angle;
}

double length_of_vector_squared(vec vec)
{
    return vec.h * vec.h + vec.k * vec.k + vec.l * vec.l;
}

double length_of_vector(vec &vec)
{
	return pow(vec.h * vec.h + vec.k * vec.k + vec.l * vec.l, 0.5);
}

vec new_vector(double h, double k, double l)
{
	vec vec;
	vec.h = h;
	vec.k = k;
	vec.l = l;

	return vec;
}

void multiply_vector(vec *vec, double mult)
{
    vec->h *= mult;
    vec->k *= mult;
    vec->l *= mult;
}

void scale_vector_to_distance(vec *vector, double new_distance)
{
	double distance = length_of_vector(*vector);

	double scale = new_distance / distance;

    multiply_vector(vector, scale);
}

double getEwaldSphereNoMatrix(vec index)
{
    if (index.l == 0)
        return 0;

    double ewald_radius = index.h * index.h + index.k * index.k
			+ index.l * index.l;

	ewald_radius /= (0 - 2 * index.l);
	double ewald_wavelength = 1 / ewald_radius;

	return ewald_wavelength;
}

double cartesian_to_distance(double x, double y)
{
    double distance = sqrt(pow(x, 2) + pow(y, 2));
    
    return distance;
}

double cartesian_to_angle(double x, double y)
{
    double angle = atan(y / x);
   
    if ((x < 0 && y > 0) ||
        (x < 0 && y < 0))
        angle += M_PI;
    
    if (x > 0 && y < 0)
        angle += M_PI * 2;
    
    return angle;
}

double cdf(double x, double mean, double sigma)
{
	double y = (x - mean) / sigma;

	double cumulative = _cdf(y);

	return cumulative;
}

double _cdf(double x)
{
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x) / sqrt(2.0);

	// A&S formula 7.1.26
	double t = 1.0 / (1.0 + p * x);
	double y = 1.0
			- (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t
					* exp(-x * x);

	return 0.5 * (1.0 + sign * y);
}

double normal_distribution(double x, double mean, double sigma)
{
	double power = 0 - pow((x - mean), 2) / (2 * sigma * sigma);
	double exp = pow(M_E, power);

	double denominator = sigma * sqrt(2 * M_PI);

	return exp / denominator;
}

double super_gaussian(double x, double mean, double sigma_0, double exponent)
{
	double correction_sigma = 2 / exponent - 1;
	double sigma = pow(M_PI / 2, correction_sigma) * sigma_0;

	double power = 0 - pow(fabs((x - mean)), exponent) / (2 * pow(sigma, exponent));

    double exp = pow(M_E, power);

	double denominator = sigma * sqrt(2 * M_PI);

	return exp / denominator;
}

double correlation_between_vectors(vector<double> *vec1,
		vector<double> *vec2)
{
	return correlation_between_vectors(vec1, vec2, NULL);
}

double correlation_between_vectors(vector<double> *vec1,
		vector<double> *vec2, vector<double> *weights)
{
	return correlation_between_vectors(vec1, vec2, weights, -1);
}

double correlation_through_origin(vector<double> *vec1,
		vector<double> *vec2, vector<double> *weights)
{
	double mean_y = weighted_mean(vec2, weights);
	double grad = gradient_between_vectors(vec1, vec2);

	double residuals_squared = 0;
	double denominator = 0;

	for (int i = 0; i < vec1->size(); i++)
	{
		double mean1 = (*vec1)[i];
		double mean2 = (*vec2)[i];

		if (mean1 != mean1 || mean2 != mean2)
			continue;

		residuals_squared += pow(mean2 - grad * mean1, 2);
		double addition = pow(mean2 - mean_y, 2);

		denominator += addition;
	}

	double R_squared = 1 - residuals_squared / denominator;

	if (R_squared < 0)
		R_squared = 0;
	double R = sqrt(R_squared);
	if (grad < 0)
		R = 0;
	if (R != R)
		R = -1;

	return R;
}

double sum(vector<double> values)
{
    double sum = 0;
    for (int i = 0; i < values.size(); i++)
    {
        sum += values[i];
    }
    
    return sum;
}

double correlation_between_vectors(vector<double> *vec1,
		vector<double> *vec2, vector<double> *weights, int exclude)
{
	double sum_x = 0;
	double sum_y = 0;
	double num = 0;
	bool needWeights = (weights != NULL);

	for (int i = 0; i < vec1->size(); i++)
	{
		if (i == exclude)
		{
			continue;
		}

		if (weights != NULL && (*weights)[i] == 0)
			continue;

		double addition = (*vec1)[i];
		if (weights != NULL)
			addition *= (*weights)[i];
		sum_x += addition;

		addition = (*vec2)[i];
		if (weights != NULL)
			addition *= (*weights)[i];
		sum_y += addition;

		double weight = (weights == NULL) ? 1 : (*weights)[i];
		num += weight;
	}
    
	double mean_x = sum_x / num;
	double mean_y = sum_y / num;

	if (mean_x != mean_x || mean_y != mean_y)
		return 0;

	double sum_x_y_minus_mean_x_y = 0;
	double sum_x_minus_mean_x_sq = 0;
	double sum_y_minus_mean_y_sq = 0;

	for (int i = 0; i < vec1->size(); i++)
	{
		if (i == exclude)
			continue;

		if (weights != NULL && (*weights)[i] == 0)
			continue;

		double addition = ((*vec1)[i] - mean_x) * ((*vec2)[i] - mean_y);
		if (weights != NULL)
			addition *= (*weights)[i];
		sum_x_y_minus_mean_x_y += addition;

		addition = std::pow((*vec1)[i] - mean_x, 2);
		if (weights != NULL)
			addition *= (*weights)[i];
		sum_x_minus_mean_x_sq += addition;

		addition = std::pow((*vec2)[i] - mean_y, 2);
		if (weights != NULL)
			addition *= (*weights)[i];

		sum_y_minus_mean_y_sq += addition;

	}

	sum_x_y_minus_mean_x_y /= num;
	sum_x_minus_mean_x_sq /= num;
	sum_y_minus_mean_y_sq /= num;

	double r = sum_x_y_minus_mean_x_y
			/ (sqrt(sum_x_minus_mean_x_sq) * sqrt(sum_y_minus_mean_y_sq));

	return r;
}

double gradient_between_vectors(vector<double> *vec1,
		vector<double> *vec2)
{

	double sum_x_y = 0;
	double sum_x_squared = 0;

	for (int i = 0; i < vec1->size(); i++)
	{
		sum_x_y += (*vec1)[i] * (*vec2)[i];
		sum_x_squared += std::pow((*vec1)[i], 2);
	}

	double grad = sum_x_y / sum_x_squared;

	return grad;
}

double least_squares_between_vectors(vector<double> *vec1,
		vector<double> *vec2, double slope)
{
	if (slope == 0)
	{
        slope = gradient_between_vectors(vec1, vec2);
	}

	if (vec1->size() == 0)
		return FLT_MAX;

	double total = 0;
	double num = 0;

	for (int i = 0; i < vec1->size(); i++)
	{
		double a = (*vec1)[i];
		double b = (*vec2)[i] * slope;

		total += pow(a - b, 2);
		num++;
	}

	return total / num;
}

double r_factor_between_vectors(vector<double> *vec1,
		vector<double> *vec2, vector<double> *weights, double scale)
{
	double sum_numerator = 0;
	double sum_denominator = 0;
	double absolute = 0;
	double allWeights = 0;

	for (int i = 0; i < vec1->size(); i++)
	{
		double weight = (*weights)[i];
		double int1 = (*vec1)[i] * scale;
		double int2 = (*vec2)[i];

		if (weight != weight || int1 != int1 || int2 != int2)
			continue;

		allWeights += weight;

		double sqrtD = 1;

		double absAddition = fabs(int1 - int2) * sqrtD;
		absAddition /= (int1 + int2) * sqrtD / 2;

		absolute += fabs(absAddition) * weight;

		sum_numerator += fabs(int1 - int2) * weight;
		sum_denominator += (int1 + int2) * weight / 2;
	}

	double r_split = sum_numerator / (sum_denominator * sqrt(2));

	return r_split;
}

double weighted_mean(vector<double> *means, vector<double> *weights)
{
	double sum = 0;
	double weight_sum = 0;

	for (int i = 0; i < means->size(); i++)
	{
		double weight = (weights == NULL ? 1 : (*weights)[i]);

		sum += (*means)[i] * weight;
		weight_sum += weight;
	}

	return sum / weight_sum;
}


bool higher(double mean1, double mean2)
{
    return mean1 > mean2;
}


double median(vector<double> *means)
{
    std::sort(means->begin(), means->end(), higher);
    size_t total = means->size();
    size_t mid = (total % 2 == 1) ? ((total - 1) / 2) : (total / 2);
    
    double midPoint = 0;
    
    if (means->size() == 0)
    {
        return std::nan(" ");
    }
    if (means->size() % 2 == 0)
    {
        midPoint = ((*means)[mid] + (*means)[mid + 1]) / 2;
    }
    else
        midPoint = (*means)[mid];
    
    
    return midPoint;
}

void histogram_gaussian(vector<double> *means, vector<int> *freq,
		double &mean, double &stdev)
{
	double sum = 0;
	double weight_sum = 0;

	for (int i = 0; i < means->size(); i++)
	{
        double weight = (*freq)[i];
		sum += (*means)[i] * weight;
		weight_sum += weight;
	}

	mean = sum / weight_sum;

	double squaredSum = 0;
	double weightSqSum = 0;

	for (int i = 0; i < means->size(); i++)
	{
        double weight = (*freq)[i];
		squaredSum += pow(mean - (*means)[i], 2) * weight;
		weightSqSum += weight;
	}

	stdev = sqrt(squaredSum / weightSqSum);
}

double standard_deviation(vector<double> *values, vector<double> *weights)
{
    double mean = weighted_mean(values, weights);
    
    return standard_deviation(values, weights, mean);
}

double standard_deviation(vector<double> *values, vector<double> *weights, double mean)
{
	double squaredSum = 0;
	double weightSqSum = 0;

	for (int i = 0; i < values->size(); i++)
	{
		double value = (*values)[i];

		if (value != value || value == FLT_MAX)
			continue;

		squaredSum += pow(mean - value, 2);
        
        double weight = (weights == NULL) ? 1 : (*weights)[i];
        
		weightSqSum += weight;
	}

	double stdev = sqrt(squaredSum / weightSqSum);

	return stdev;
}

void regression_line(vector<boost::tuple<double, double, double> > values, double &intercept, double &gradient)
{
	double sigma_x = 0;
	double sigma_y = 0;
	double sigma_x_y = 0;
	double sigma_x_2 = 0;
	double weight_n = 0;

	for (int i=0; i < values.size(); i++)
	{
		double x = boost::get<0>(values[i]);
		double y = boost::get<1>(values[i]);
		double weight = boost::get<2>(values[i]);

		sigma_x += x * weight;
		sigma_y += y * weight;
		sigma_x_y += x * y * weight;
		sigma_x_2 += x * x * weight;
		weight_n += weight;
	}

	double mean_x = sigma_x / weight_n;
	double mean_y = sigma_y / weight_n;

	double sxy = sigma_x_y - sigma_x * sigma_y / weight_n;
	double sxx = sigma_x_2 - pow(sigma_x, 2) / weight_n;

	gradient = sxy / sxx;

	intercept = mean_y - gradient * mean_x;
}

double minimizeParam(double &step, double &param, double (*score)(void *object),
                     void *object)
{
    return minimizeParameter(step, &param, score, object);
}

double minimizeParameter(double &step, double *param, double (*score)(void *object),
                                    void *object)
{
    double param_trials[3];
    double param_scores[3];
    
    int j = 0;
    int param_min_num = 1;
    
    double bestParam = *param;
    
    for (double i = bestParam - step; j < 3; i += step)
    {
        *param = i;
        
        double aScore = (*score)(object);
        
        if (aScore != aScore) aScore = FLT_MAX;
        
        param_scores[j] = aScore;
        param_trials[j] = i;
        j++;
    }
    
    double param_min_score = param_scores[1];
    
    for (int i = 0; i < 3; i++)
        if (param_scores[i] < param_min_score)
        {
            param_min_score = param_scores[i];
            param_min_num = i;
        }
    
    *param = param_trials[param_min_num];
    (*score)(object);
    
    if (param_min_num == 1)
        step /= 2;
    
    return param_scores[param_min_num];
}

std::map<double, int> histogram(std::vector<double> values, double step)
{
    std::map<int, int> tempHistogram;
    
    for (int i = 0; i < values.size(); i++)
    {
        if (values[i] != values[i])
            continue;
        
        int category = values[i] / step;
        
        if (values[i] == 0)
            category = 0;
        
        tempHistogram[category]++;
    }
    
    std::map<double, int> realHistogram;
    
    if (values.size() == 0)
        return realHistogram;
    
   // Logger::mainLogger->addString("Histogram has " + std::to_string((int)tempHistogram.size()) + " categories.");
    
    int minCategory = INT_MAX;
    int maxCategory = -INT_MAX;
    
    for (std::map<int, int>::iterator it = tempHistogram.begin(); it != tempHistogram.end(); it++)
    {
        int category = it->first;
        
        if (category < minCategory)
            minCategory = category;
        if (category > maxCategory)
            maxCategory = category;
    }
    
    for (int i = minCategory; i <= maxCategory; i++)
    {
        double realStep = i * step;
        int frequency = tempHistogram[i];
        
        realHistogram[realStep] = frequency;
    }
    
    return realHistogram;
}

void histogramCSV(std::string filename, std::map<double, int> map1, std::map<double, int> map2)
{
    std::ofstream stream;
    stream.open(filename.c_str());
    
    for (std::map<double, int>::iterator it = map1.begin(); it != map1.end(); it++)
    {
        double category = it->first;
        int first = it->second;
        int second = 0;
        
        if (map2.count(category))
        {
            second = map2[category];
        }
        
        stream << category << "," << first << "," << second << std::endl;
    }
    
    stream.close();
}

void printDesc(vec hkl)
{
    Logger::mainLogger->addString(desc(hkl));
}

std::string prettyDesc(vec hkl)
{
    std::ostringstream logged;
    logged << "(" << hkl.h << ", " << hkl.k << ", " << hkl.l << ")";
    return logged.str();
}

std::string desc(vec hkl)
{
    std::ostringstream logged;
    logged << hkl.h << "\t" << hkl.k << "\t" << hkl.l;
    return logged.str();
}