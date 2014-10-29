#include <math.h>

void polygon_add(float* polygon, float* vect, int q_vertices, float* ret)
{
	unsigned i;
	for (i=0; i<q_vertices*2; i+=2)
	{
		ret[i]   = polygon[i]   + vect[0];
		ret[i+1] = polygon[i+1] + vect[1];
	}
}


void polygon_sub(float* polygon, float* vect, int q_vertices, float* ret)
{
	unsigned i;
	for (i=0; i<q_vertices*2; i+=2)
	{
		ret[i]   = polygon[i]   - vect[0];
		ret[i+1] = polygon[i+1] - vect[1];
	}
}

void polygon_rotate(float* polygon, float angle, float* origin, int q_vertices, float* ret)
{
	unsigned i;
	for (i=0; i<q_vertices*2; i+=2)
	{
		float x = polygon[i];
		float y = polygon[i+1];
		float ox = origin[0];
		float oy = origin[1];
		float ca = cos(angle);
		float sa = sin(angle);
		ret[i]   = ((x - ox) * ca) - ((y - oy) * sa);
		ret[i+1] = ((x - ox) * sa) + ((y - oy) * ca);
	}
}
