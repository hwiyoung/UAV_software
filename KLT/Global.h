/********************************************************
 * Project:			Automated AT
 * Last updated:	17 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _GLOBAL_H_
#define	_GLOBAL_H_

#pragma once

#ifndef	M_PI
#define	M_PI 3.14159265358979323846
#endif

#define	DEG_H_CIRCLE	180
#define	DEG_TO_RAD		(M_PI / DEG_H_CIRCLE)
#define	RAD_TO_DEG		(DEG_H_CIRCLE / M_PI)

double deg2rad(double degrees);

double rad2deg(double radians);

#endif // _GLOBAL_H_
