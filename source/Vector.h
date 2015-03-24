/*
 * Vector.h
 *
 *  Created on: 27 Aug 2014
 *      Author: helenginn
 */

#ifndef VECTOR_H_
#define VECTOR_H_
#include <vector>

#include <boost/tuple/tuple.hpp>

typedef struct
{
	double h;
	double k;
	double l;
} vec;

vec cross_product_for_vectors(vec vec1, vec vec2);
vec new_vector(double h, double k, double l);
double length_of_vector(vec vect);
double angleBetweenVectors(vec vec1, vec vec2);
vec copy_vector(vec old_vec);
void add_vector_to_vector(vec *vec1, vec vec2);
vec vector_between_vectors(vec vec1, vec vec2);
void take_vector_away_from_vector(vec vec1, vec *vec2);
void scale_vector_to_distance(vec *vec, double new_distance);
double getEwaldSphereNoMatrix(vec index);
double cdf(double x, double mean, double sigma);
double _cdf(double x);
double normal_distribution(double x, double mean, double sigma);
double super_gaussian(double x, double mean, double sigma_0, double exponent);

double minimizeParam(double &step, double &param, double (*score)(void *object),
                         void *object);
double minimizeParameter(double &step, double &param, double (*score)(void *object),
                       void *object);

void regression_line(std::vector<boost::tuple<double, double, double> > values, double &intercept, double &gradient);
double correlation_between_vectors(std::vector<double> *vec1,
		std::vector<double> *vec2, std::vector<double> *weights, int exclude);
double correlation_between_vectors(std::vector<double> *vec1,
		std::vector<double> *vec2, std::vector<double> *weights);
double correlation_between_vectors(std::vector<double> *vec1,
		std::vector<double> *vec2);
double correlation_through_origin(std::vector<double> *vec1,
		std::vector<double> *vec2, std::vector<double> *weights = NULL);
double least_squares_between_vectors(std::vector<double> *vec1,
		std::vector<double> *vec2, double slope);
double gradient_between_vectors(std::vector<double> *vec1,
		std::vector<double> *vec2);
double minimize_gradient_between_vectors(std::vector<double> *vec1,
		std::vector<double> *vec2);
double weighted_mean(std::vector<double> *means, std::vector<double> *weights = NULL);
void histogram_gaussian(std::vector<double> *means, std::vector<int> *freq, double &mean, double &stdev);
double least_squares_gaussian_fit(std::vector<double> *means,
		std::vector<int> *freq);
double standard_deviation(std::vector<double> *values, std::vector<double> *weights = NULL);
double r_factor_between_vectors(std::vector<double> *vec1,
		std::vector<double> *vec2, std::vector<double> *weights, double scale);

double cartesian_to_distance(double x, double y);
double cartesian_to_angle(double x, double y);

#endif /* VECTOR_H_ */
