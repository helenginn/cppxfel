/*
 * Index.h
 *
 *  Created on: 26 Dec 2014
 *      Author: helenginn
 */

#ifndef INDEX_H_
#define INDEX_H_

#include "parameters.h"

class Index
{
public:
	Index();
	virtual ~Index();
};

void index(char **argv, int argc);
ShoeboxPtr newShoebox();

#endif /* INDEX_H_ */
